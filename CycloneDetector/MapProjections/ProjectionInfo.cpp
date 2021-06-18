#include "./ProjectionInfo.h"

#include <limits>
#include <algorithm>
#include <cstring>
#include <cassert>

#include "./Projections/Mercator.h"
#include "./Projections/LambertConic.h"
#include "./Projections/Equirectangular.h"

#include "MapProjectionUtils.h"

using namespace Projections;

/// <summary>
/// Get name of projection
/// </summary>
/// <returns></returns>
template <typename Proj>
const char* ProjectionInfo<Proj>::GetName() const
{
	return static_cast<const Proj*>(this)->GetNameInternal();
}

//=======================================================================
// Main interface
//=======================================================================

template <typename Proj>
ProjectionInfo<Proj>::ProjectionInfo(PROJECTION curProjection)
	: IProjectionInfo(curProjection)		
{	

}

template <typename Proj>
std::tuple<double, double, double, double> ProjectionInfo<Proj>::GetFrameBotLeftTopRight(
	const Coordinate & botLeft, const Coordinate & topRight)
{
	ProjectedValue tmpMinPixel = static_cast<Proj*>(this)->ProjectInternal(botLeft);
	ProjectedValue tmpMaxPixel = static_cast<Proj*>(this)->ProjectInternal(topRight);

	ProjectedValue minVal, maxVal;
	
	minVal.x = std::min(tmpMinPixel.x, tmpMaxPixel.x);
	minVal.y = std::min(tmpMinPixel.y, tmpMaxPixel.y);

	maxVal.x = std::max(tmpMinPixel.x, tmpMaxPixel.x);
	maxVal.y = std::max(tmpMinPixel.y, tmpMaxPixel.y);

	return std::make_tuple(minVal.x, minVal.y, maxVal.x, maxVal.y);
}


template <typename Proj>
void ProjectionInfo<Proj>::SetFrame(const ProjectionFrame & frame)
{	
	this->frame.h = frame.h;
	this->frame.w = frame.w;
	this->frame.hAR = frame.hAR;
	this->frame.wAR = frame.wAR;
	this->frame.hPadding = frame.hPadding;
	this->frame.wPadding = frame.wPadding;
	
	this->frame.projPrecomX = frame.projPrecomX;
	this->frame.projPrecomY = frame.projPrecomY;
	this->frame.stepType = frame.stepType;

	this->frame.repeatNegCount = frame.repeatNegCount;
	this->frame.repeatPosCount = frame.repeatPosCount;
	this->frame.ww = frame.ww;
	this->frame.hh = frame.hh;

	//set some default values from old frame,
	//but they will be recomputed in ComputeAABB
	this->frame.min = frame.min;
	this->frame.max = frame.max;

	
	this->ComputeAABB(this->frame.min, this->frame.max);
}



/// <summary>
/// Set current data active frame based on image bottom left and
/// top right coordinate
/// 
/// This assumes that data are plotted in 2D image and image has corners
/// These corners do not have to correspond to AABB of coordinates
/// For example: Lambert - image is square but corners are not AABB
/// 
/// If width or heigth is zero, the size is determined based on
/// other dimension. In this case, it is always expected to keepAR
/// 
/// StepType - determine, where coordinate botLeft and topRight is positioned within the pixel
/// It can be at pixel center or on its border
/// First position [0, 0] is always correct,
/// Max position [w, h] can be shifted by -1 based on mode:
/// 
/// If we have data:
/// Border mode: 
/// width = 3 pixels (a---b---c---d), GPS is located at a, b, c, d
/// In this case, there are 3 "steps" between GPS positions, but 4 positions in total
/// For this, max position is calculated directly from width that correspond to steps
/// a + 3 * '---' => d
/// 
/// Center mode:
/// width = 3 pixels (|-a-|-b-|-c-|), GPS is located at a, b, c
/// In this case, there are only 2 steps between GPS positions, but 3 positions in total
/// For this, max position is calculated by subtract 1 from width to get correct last.
/// a + (3 - 1) * '-x-' => c
/// 
/// </summary>
/// <param name="botLeft">Bottom left coordinate</param>
/// <param name="topRight">Top right coordinate</param>
/// <param name="w">frame width</param>
/// <param name="h">frame height</param>
/// <param name="stepType">type of pixel steps</param> 
/// <param name="keepAR">keep AR of data (default: true) 
/// if yes, data are enlarged beyond AABB to keep AR
/// </param>
template <typename Proj>
void ProjectionInfo<Proj>::SetRawFrame(const Coordinate & botLeft, const Coordinate & topRight,
	MyRealType w, MyRealType h, STEP_TYPE stepType, bool keepAR)
{		
	//calculate minimum internal projection value				

	auto[minX, minY, maxX, maxY] = static_cast<Proj*>(this)->GetFrameBotLeftTopRight(botLeft, topRight);

	this->frame.stepType = stepType;

	//Calculate width / height of internal projection		
	MyRealType projW = maxX - minX;
	MyRealType projH = maxY - minY;

	//----------------------------------------------------------

	//calculate scale in width and height		
	if (w != 0)
	{
		this->frame.w = w;
	}
	else
	{
		MyRealType wAr = projW / projH;
		this->frame.w = h * wAr;
	}

	if (h != 0)
	{
		this->frame.h = h;		
	}
	else 
	{
		MyRealType hAr = projH / projW;
		this->frame.h = w * hAr;
	}

	this->frame.wAR = (this->frame.w - stepType) / projW;
	this->frame.wPadding = 0;

	this->frame.hAR = (this->frame.h - stepType) / projH;
	this->frame.hPadding = 0;

	// Using different ratios for width and height will cause the map to be stretched,
	// but fitting the desired region
	// If we want to keep AR of the map, the AR will be correct, but frame will no
	// be the one, we have set. One dimmension will be changed to keep AR	
	if (keepAR)
	{
		MyRealType globalAR = std::min(this->frame.wAR, this->frame.hAR);
		this->frame.wAR = globalAR;
		this->frame.hAR = globalAR;

		this->frame.wPadding = ((this->frame.w - stepType) - (this->frame.wAR * projW)) * MyRealType(0.5);
		this->frame.hPadding = ((this->frame.h - stepType) - (this->frame.hAR * projH)) * MyRealType(0.5);
	}

	this->frame.projPrecomX = -this->frame.wPadding + this->frame.wAR * minX;
	this->frame.projPrecomY = -(this->frame.h - stepType) + this->frame.hPadding - this->frame.hAR * minY;

	this->frame.min = botLeft;
	this->frame.max = topRight;

	this->CalculateWrapRepeat(botLeft, topRight);
}

/// <summary>
/// Calculate wrap around counts
/// - How many times are coordinates wrapped around the world
/// in longitude
/// - for orthogonal projections
/// </summary>
/// <param name="botLeft"></param>
/// <param name="topRight"></param>
template <typename Proj>
void ProjectionInfo<Proj>::CalculateWrapRepeat(const Coordinate& botLeft, const Coordinate& topRight)
{
	if (Proj::ORTHOGONAL_LAT_LON == false)
	{
		//cannot repeat for non-orthogonal lat/lon projections
		return;
	}

	//calculate repeat count
	this->frame.repeatNegCount = 0;
	this->frame.repeatPosCount = 0;

	MyRealType lon = botLeft.lon.deg();
	if (lon < -180.0)
	{
		while (lon < -360.0)
		{
			this->frame.repeatNegCount++;
			lon += 360.0;
		}

		//map range to 0 - 1
		lon += 360.0;
		MyRealType in01 = (lon - -180.0) * (1.0) / 360.0;

		this->frame.repeatNegCount += (1.0 - in01);
	}

	lon = topRight.lon.deg();
	if (lon > 180.0)
	{
		while (lon > 360.0)
		{
			this->frame.repeatPosCount++;
			lon -= 360.0;
		}

		//map range to 0 - 1
		lon -= 360.0;
		MyRealType in01 = (lon - -180.0) * (1.0) / 360.0;

		this->frame.repeatPosCount += in01;
	}

	Projections::Coordinate bbMin, bbMax;

	bbMin.lat = -90.0_deg;
	bbMin.lon = -180.0_deg;
	
	bbMax.lat = 90.06_deg;
	bbMax.lon = 180.0_deg;
	
	Pixel<MyRealType> pp1 = this->Project<MyRealType>(bbMin);
	Pixel<MyRealType> pp2 = this->Project<MyRealType>(bbMax);

	this->frame.ww = pp2.x - pp1.x;
	this->frame.hh = pp1.y - pp2.y;

}

/// <summary>
/// Same logic as SetRawFrame
/// BUT: min/max coordinate of frame is calculated as AABB 
/// using inverse projection
/// </summary>
/// <param name="botLeft"></param>
/// <param name="topRight"></param>
/// <param name="w"></param>
/// <param name="h"></param>
/// <param name="stepType"></param>
/// <param name="keepAR"></param>
template <typename Proj>
void ProjectionInfo<Proj>::SetFrame(const Coordinate & botLeft, const Coordinate & topRight,
	MyRealType w, MyRealType h, STEP_TYPE stepType, bool keepAR)
{				
	this->SetRawFrame(botLeft, topRight, w, h, stepType, keepAR);
	this->ComputeAABB(this->frame.min, this->frame.max);
}


/// <summary>
/// Set current data active frame based on AABB
/// For projections that have orthogonal lat/lon,
/// this is same as SetFrame, since botLeft = min and topRight = max
/// 
/// However, for other projections like Lambert, this will behave differently
/// 
/// </summary>
/// <param name="min"></param>
/// <param name="max"></param>
/// <param name="w"></param>
/// <param name="h"></param>
/// <param name="keepAR"></param>
template <typename Proj>
void ProjectionInfo<Proj>::SetFrameFromAABB(const Coordinate & min, const Coordinate & max,
	MyRealType w, MyRealType h, STEP_TYPE stepType, bool keepAR)
{
	//not working correctly for non-orthogonal lat/lon projections !!!!!
	//todo
	this->SetFrame(min, max, w, h, stepType, keepAR);
}

/// <summary>
/// calculate how many "degrees" is for one "pixel" in lat / lon
/// cordinates are boundaries of pixels :
/// This is correct only for projections that have orthogonal lat / lon
/// (like Mercator), but not for eg. Lambert
/// 
/// Eg:
/// For 
/// Image with frame 2x2 px and border mode:
///    -180  +180
///  +90 *-----*
///      |  |  |
///      -------
///      |  |  |
///  -90 *-----*
/// Step latitude : (90 - (-90)) / 2 = 90
/// Step longitude : (180 - (-180)) / 2 = 180
/// </summary>
/// <returns></returns>
template <typename Proj>
Coordinate ProjectionInfo<Proj>::GetDeltaStep() const
{	
	Coordinate step;
	step.lat = Latitude::rad((this->frame.max.lat.rad() - this->frame.min.lat.rad()) / 
		(this->frame.h - this->frame.stepType));
	step.lon = Longitude::rad((this->frame.max.lon.rad() - this->frame.min.lon.rad()) / 
		(this->frame.w - this->frame.stepType));

	return step;
}

/// <summary>
/// Get projection top left corner
/// </summary>
/// <returns></returns>
template <typename Proj>
Coordinate ProjectionInfo<Proj>::GetTopLeftCorner() const
{
	return this->ProjectInverse({ 0, 0 });
}

template <typename Proj>
const ProjectionFrame & ProjectionInfo<Proj>::GetFrame() const
{
	return this->frame;
}


/// <summary>
/// Get "pixels" on line. For each "pixel", callback function is called
/// </summary>
/// <param name="start"></param>
/// <param name="end"></param>
/// <param name="callback"></param>
template <typename Proj>
void ProjectionInfo<Proj>::LineBresenham(Pixel<int> start, Pixel<int> end,
	std::function<void(int x, int y)> callback) const
{
	if ((start.x < 0) || (start.y < 0))
	{
		return;
	}

	if ((end.x < 0) || (end.y < 0))
	{
		return;
	}

	
	//for coordinates at pixel corner:
	//[0, 0] is at pixel corner
	//[w, h] is at pixel corner
	//for coordinates at pixel center:
	//use [w - 1, h - 1] (it correspond to [0, 0] for pixel center)
	int ww = static_cast<int>(this->frame.w) - this->frame.stepType;
	int hh = static_cast<int>(this->frame.h) - this->frame.stepType;


	if ((start.x > ww) || (start.y > hh))
	{
		return;
	}

	if ((end.x > ww) || (end.y > hh))
	{
		return;
	}
	



	int dx = std::abs(end.x - start.x);
	int dy = std::abs(end.y - start.y);
	int sx, sy, e2;


	(start.x < end.x) ? sx = 1 : sx = -1;
	(start.y < end.y) ? sy = 1 : sy = -1;
	int err = dx - dy;

	while (1)
	{		
		callback(start.x, start.y);
		if ((start.x == end.x) && (start.y == end.y))
		{
			break;
		}
		e2 = 2 * err;
		if (e2 > -dy)
		{
			err = err - dy;
			start.x = start.x + sx;
		}
		if (e2 < dx)
		{
			err = err + dx;
			start.y = start.y + sy;
		}
	}
}



/// <summary>
/// Compute AABB for current active frame
/// It uses Bresenham lines around the image 
/// E.g.: [0,0] -> [0, h] and reproject each pixel, and from those, AABB is calculated
/// </summary>
/// <param name="min"></param>
/// <param name="max"></param>
template <typename Proj>
void ProjectionInfo<Proj>::ComputeAABB(Coordinate & min, Coordinate & max) const
{	
	//for coordinates at pixel corner:
	//[0, 0] is at pixel corner
	//[w, h] is at pixel corner
	//for coordinates at pixel center:
	//use [w - 1, h - 1] (it correspond to [0, 0] for pixel center)
	int ww = static_cast<int>(this->frame.w) - this->frame.stepType;
	int hh = static_cast<int>(this->frame.h) - this->frame.stepType;

			
	std::vector<Coordinate> border;

	if (static_cast<const Proj*>(this)->ORTHOGONAL_LAT_LON)
	{
		Coordinate c = this->ProjectInverse({ 0, 0 });
		border.push_back(std::move(c));

		c = this->ProjectInverse({ ww, hh });
		border.push_back(std::move(c));
	}
	else
	{
		this->LineBresenham({ 0,0 }, { 0, hh },
			[&](int x, int y) -> void {
				Coordinate c = this->ProjectInverse({ x, y });
				border.push_back(std::move(c));
			});
		this->LineBresenham({ 0,0 }, { ww, 0 },
			[&](int x, int y) -> void {
				Coordinate c = this->ProjectInverse({ x, y });
				border.push_back(std::move(c));
			});
		this->LineBresenham({ ww, hh }, { 0, hh },
			[&](int x, int y) -> void {
				Coordinate c = this->ProjectInverse({ x, y });
				border.push_back(std::move(c));
			});
		this->LineBresenham({ ww, hh }, { 0, hh },
			[&](int x, int y) -> void {
				Coordinate c = this->ProjectInverse({ x, y });
				border.push_back(std::move(c));
			});
	}
	
	ProjectionUtils::ComputeAABB(border, min, max);
}




//=====

template class Projections::ProjectionInfo<LambertConic>;
template class Projections::ProjectionInfo<Mercator>;
template class Projections::ProjectionInfo<Equirectangular>;
