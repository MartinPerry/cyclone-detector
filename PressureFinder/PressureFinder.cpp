#include "./PressureFinder.h"

#include <stack>
#include <array>
#include <unordered_set>

#include "./ADT/BVH.h"

#include "./RasterData/ImageGausFilter.h"
#include "./RasterData/ImageUtils.h"
#include "./RasterData/MarchingSquares.h"

#include "./Graphics/2d/PolyLine.h"
#include "./Graphics/2d/Polygon.h"

#include "./Math/MathUtils.h"
#include "./Math/Vector2.h"

#include "./Utils/Logger.h"

#include "./MapProjections/IProjectionInfo.h"
#include "./MapProjections/ProjectionInfo.h"
#include "./MapProjections/MapProjectionUtils.h"
#include "./MapProjections/Projections/Equirectangular.h"
#include "./MapProjections/Projections/Mercator.h"
#include "./MapProjections/Projections/LambertConic.h"

#include "PressureExtrema.h"

std::string PressureFinder::debugName = "";
std::string PressureFinder::debugDirName = "";


/// <summary>
/// ctor
/// Load data from RAW file with fileName
/// (RAW file stores directly float values)
/// </summary>
/// <param name="fileName"></param>
/// <param name="w"></param>
/// <param name="h"></param>
PressureFinder::PressureFinder(const char * fileName, int w, int h) :
	//rawData(Image2d<float>::CreateFromRawFile(1097, 657, ColorSpace::PixelFormat::GRAY, fileName)),
	rawData(Image2d<float>::CreateFromRawFile(w, h, ColorSpace::PixelFormat::GRAY, fileName)),
	latCorrectionFactor(15),
	maskRadius(20),
	contourStep(200),
	minCentersDistanceThresholdPixel(22),
	smallAreaThreshold(AreaThreshold::CreatePixels(100)),
	extremaCorrectRatio(0.6)
{
	rawData.FindMinMax(0, minValue, maxValue);
	//rawData.RunGauss(3);
}

/// <summary>
/// ctor
/// </summary>
/// <param name="data"></param>
/// <param name="w"></param>
/// <param name="h"></param>
PressureFinder::PressureFinder(const std::vector<float> & data, int w, int h) :
	rawData(Image2d<float>(w, h, data, ColorSpace::PixelFormat::GRAY)),	
	latCorrectionFactor(15),
	maskRadius(20),
	contourStep(200),
	minCentersDistanceThresholdPixel(22),
	smallAreaThreshold(AreaThreshold::CreatePixels(100)),
	extremaCorrectRatio(0.6)
{
	rawData.FindMinMax(0, minValue, maxValue);		
}

/// <summary>
/// ctor
/// </summary>
/// <param name="input"></param>
PressureFinder::PressureFinder(Image2d<float>&& input) :
	rawData(std::move(input)),
	latCorrectionFactor(15),
	maskRadius(20),
	contourStep(200),
	minCentersDistanceThresholdPixel(22),
	smallAreaThreshold(AreaThreshold::CreatePixels(100)),
	extremaCorrectRatio(0.6)
{
	rawData.FindMinMax(0, minValue, maxValue);
}


/// <summary>
/// dtor
/// </summary>
PressureFinder::~PressureFinder()
{
}

double PressureFinder::GetContoursStep() const
{
	return this->contourStep;
}

int PressureFinder::GetSystemsCount() const
{
	int count = 0;
	for (auto& c : this->extremas)
	{
		if (c.IsValid() == false)
		{
			continue;
		}
		count++;
	}

	return count;
}

/// <summary>
/// Get all detected extremas
/// </summary>
/// <returns></returns>
const std::vector<PressureExtrema>& PressureFinder::GetExtremas() const
{
	return this->extremas;
}

const std::vector<Contour>& PressureFinder::GetContours() const
{
	return this->contours;
}

const Image2d<float>& PressureFinder::GetRawData() const
{
	return this->rawData;
}

/// <summary>
/// Get minimal value from data
/// </summary>
/// <returns></returns>
float PressureFinder::GetDataMin() const
{
	return this->minValue;
}

/// <summary>
/// Get maximal value from data
/// </summary>
/// <returns></returns>
float PressureFinder::GetDataMax() const
{
	return this->maxValue;
}

/// <summary>
/// Correction factor for latitude when areas size is calculated
/// On equator, area is taken directly from smallAreaThreshold
/// On poles, it is multipled by this factor
/// and between, the factor is linearly interpolated
/// 
/// Default: 15
/// 
/// To disable this feature, set factor to 0.0
/// </summary>
/// <param name="f"></param>
void PressureFinder::SetAreaLatitudeCorrectionFactor(double f)
{
	this->latCorrectionFactor = f;
}

/// <summary>
/// For each extrema center find the closest center in another
/// extrema. If centers are too close, remove the center with smaller area
/// 
/// Default: 500
/// 
/// To disable this feature, set value to 0
/// </summary>
/// <param name="v"></param>
void PressureFinder::SetMinCentersDistanceThresholdPixel(int v)
{
	this->minCentersDistanceThresholdPixel = v;
}

/// <summary>
/// Area on equator aprox. 20'000 km^2 = 10 * 10 px for ICON 2880x1441
/// as we go to the poles, area with the same pixel size decreases
/// in equirectangular projection
/// 
/// Default: 10*10px
/// </summary>
/// <param name="v"></param>
void PressureFinder::SetSmallAreaThreshold(AreaThreshold v)
{
	this->smallAreaThreshold = v;
}

/// <summary>
/// pixel is considered extrema if ratio of neighborhood pixels
/// that correspond to mask against "incorrect" pixels is above extremaCorrectRatio
/// (i.e. extremaCorrectRatio of area pixels is correct)
/// 
/// Default: 0.6
/// </summary>
/// <param name="v"></param>
void PressureFinder::SetExtremaCorrectRatio(double v)
{
	this->extremaCorrectRatio = v;
}

/// <summary>
/// Set mask radius
/// Mask is used to detect extrema
/// 
/// If set to 0 = mask radius is auto calculated
/// </summary>
/// <param name="r"></param>
void PressureFinder::SetExtremaMaskRadius(int r)
{
	this->maskRadius = r;	
}

/// <summary>
/// Set step between contours in Pa
/// </summary>
/// <param name="step"></param>
void PressureFinder::SetContoursStep(double step)
{
	this->contourStep = step;
}

/// <summary>
/// Run the algorithm
/// </summary>
void PressureFinder::Run()
{
	this->contours.clear();
	this->extremas.clear();

	MY_LOG_INFO("Data min: %f, Data max: %f", minValue, maxValue);

	this->FindContours();

	this->FindExtrema();

#ifndef DISABLE_DEBUG_OUTPUT
	this->RenderContours();
#endif

	
	int count = this->GetSystemsCount();	
	MY_LOG_INFO("Detected systems: %d", count);
}

/// <summary>
/// Find contours in image
/// using Marching Squares
/// </summary>
void PressureFinder::FindContours()
{
	//build contours
	for (double t = minValue; t <= maxValue; t += this->contourStep)
	{
		MarchingSquares ms;

		//set "border" value for rawdata
		//we use the condition to enclose contours that go outside the image
		//and should be wrapped around
		ms.SetBorderValue((t < 0.75 * maxValue) ? maxValue : minValue);
		
		ms.Run<float>(rawData, static_cast<float>(t));

		for (const auto & l : ms.GetAllContours())
		{
			if (l.Size() < 10)
			{
				//contour has too few points (in this case point = pixel)
				continue;
			}

			//auto lSimple = l.CreateSimplified(5);			
			//contours.push_back(d2::Polygon::CreateFromPolyline(lSimple));

			//to speed up further calculations, simplify contour polyline			
			contours.emplace_back(t, d2::Polygon::CreateFromPolyline(l.CreateSimplified(10)));
		}

	}

	this->BuildContoursHierarchy();
}

/// <summary>
/// Build contours hierarchy
/// - find extrema contours (no contour is inside)
/// - find contour parents
/// - find contour childrens
/// - remove small extrema contours and "move" the extrema
/// to its parent
/// </summary>
void PressureFinder::BuildContoursHierarchy()
{
	//build extrema contours
	for (size_t i = 0; i < contours.size(); i++)
	{
		Contour & c = contours[i];
			
		for (size_t j = 0; j < contours.size(); j++)
		{
			if (i == j) continue;
			
			//test if contour has some other points inside
			//if not - extrema
			if (c.c.IsPointInside(contours[j].c[0]))
			{								
			
				//build childrens
				//contour[j] is inside c => it is its children
				c.childCount++;				
				c.childs.push_back(&contours[j]);
			}
		}		

		if (c.childCount == 0)
		{				
			c.isExtrema = true;			
		}
	}

#ifndef DISABLE_DEBUG_OUTPUT
	RenderContours(contours);
#endif

	//build parents		
	for (size_t i = 0; i < contours.size(); i++)
	{
		Contour & c = contours[i];

		int minParentChildCount = std::numeric_limits<int>::max();
		Contour * parent = nullptr;
		
		for (size_t j = 0; j < contours.size(); j++)
		{
			if (i == j) continue;
						
			//test ic contour c is inside contour[j] 
			if (contours[j].c.IsPointInside(c.c[0]))
			{
				//contour is inside => contour[j] is one if its "parent"
				//from all "parents", the direct parent has minimal number of childrens
				//because each outter parent has at least +1 children
				
				if (minParentChildCount > contours[j].childCount)
				{
					parent = &contours[j];
					minParentChildCount = parent->childCount;
				}
			}
		}

		c.parent = parent;
	}
	
	//filter extremas
		
	//number of childrens depends on step between contours !!!!
	//we use inner area of 30 hPa as a threshold
	//(the quite high threshold is used to correctly process small areas with a large pressure change)
	//const int CHILDS_COUNT_THRESHOLD = static_cast<int>(3000 / this->contourStep);
	
	bool isChanged = false;
	do
	{
		isChanged = false;

		for (size_t i = 0; i < contours.size(); i++)
		{

			Contour & c = contours[i];

			if (c.isExtrema == false)
			{
				continue;
			}
			
			double area = 0;
			double areaThreshold = 0;

			if (this->smallAreaThreshold.unit == AreaThreshold::Unit::Km2)
			{
				area = this->CalcArea(c, this->smallAreaThreshold.proj);
				areaThreshold = this->smallAreaThreshold.value;
			}
			else
			{
				const auto& bb = c.c.GetBoundingBox();
				areaThreshold = this->GetLatitudeAreaCorrection(bb.GetCentroid().y);
				area = bb.GetArea();
			}

			//check if contour has parrent and contour area is small enough
			if ((c.parent) && (area < areaThreshold))
			{
				isChanged = true;
				c.isExtrema = false; //current contour is not extrema (too small)

				
				double areaParent = 0;
				if (this->smallAreaThreshold.unit == AreaThreshold::Unit::Km2)
				{
					areaParent = this->CalcArea(*c.parent, this->smallAreaThreshold.proj);
				}
				else
				{
					areaParent = c.parent->c.GetBoundingBox().GetArea(); //get parent area
				}
				

				//area if parent is not "large"
				//it can be up to 100x larger than child threshold
				if (areaParent < areaThreshold * 100)
				{
					//check if parent has any children extrema
					//different than current contour c
					bool hasOtherChildExtrema = false;
					for (Contour * child : c.parent->childs)
					{
						if ((child != &c) && (child->isExtrema))
						{
							//there is another extrema
							hasOtherChildExtrema = true;
							break;
						}
					}

					if (hasOtherChildExtrema == false)
					{
						c.parent->isExtrema = true;
					}
				}

				//we have to restart main loop and iterate contrours again
				//because now parent can be set as an extreme
				//and parent can be on lower index i
				//now we apply same condition on parent extrema
				//and this may cause parent extrema to enlarge egain to its parent
				
				//NOTE: without this, it should work mostly as well since
				//parent is taken only if "large" enough
				break;
			}
		}
	} while (isChanged);
	


	double totalSquareRootArea = 0;
	
	//store only extremas in special array
	for (size_t i = 0; i < contours.size(); i++)
	{
		Contour & c = contours[i];

		if (c.isExtrema)
		{
			totalSquareRootArea += sqrt(c.c.CalcArea());			

			extremas.emplace_back(c);
		}
	}	


	if (this->maskRadius == 0)
	{
		//no mask radius set
		//use auto-calculated raius
		this->SetExtremaMaskRadius(static_cast<int>(totalSquareRootArea / extremas.size()));

		if (this->maskRadius < 10)
		{
			//mask radius is too small
			this->maskRadius = 10;
		}

		MY_LOG_INFO("Autocalculated mask radius: %d\n", this->maskRadius);		
	}
	
}

/// <summary>
/// Calculate contour area
/// - if proj is not set or is uknown - return value in pixels
/// - if projection is set, calculate area of contour AABB in km^2
/// </summary>
/// <param name="c"></param>
/// <param name="proj"></param>
/// <returns></returns>
double PressureFinder::CalcArea(const Contour& c, Projections::IProjectionInfo* proj)
{	
	if (proj == nullptr)
	{
		return c.c.CalcArea();
	}

	if (auto eq = dynamic_cast<Projections::Equirectangular*>(proj))
	{
		return CalcAreaProj(c, eq);
	}
	if (auto me = dynamic_cast<Projections::Mercator*>(proj))
	{
		return CalcAreaProj(c, me);
	}
	if (auto lc = dynamic_cast<Projections::LambertConic*>(proj))
	{
		return CalcAreaProj(c, lc);
	}
	
	return c.c.CalcArea();
}

template <typename Proj>
double PressureFinder::CalcAreaProj(const Contour& c, Proj* proj)
{
	const auto& aabb = c.c.GetBoundingBox();

	std::vector<Projections::Pixel<float>> px;
	for (int i = 0; i < 4; i++)
	{
		auto co = aabb.GetCorner(i);
		px.push_back({ float(co.x), float(co.y) });
	}

	for (size_t i = 0; i < c.c.Size(); i++)
	{
		//px.push_back({ float(c.c[i].x), float(c.c[i].y) });
	}

	return Projections::ProjectionUtils::CalcArea(px, proj) / (1000.0 * 1000.0);
}

/// <summary>
/// calculate latitude area correction
/// on equator, we have base area given by smallAreaThreshold
/// as we go near the pole, the same pixel size area correspond to a smaller areas
/// we need to keep area more consistent, so we use correction factor
/// </summary>
/// <param name="y"></param>
/// <returns></returns>
double PressureFinder::GetLatitudeAreaCorrection(double y)
{		
	if (latCorrectionFactor == 0)
	{
		//correczion factor is disabled
		return smallAreaThreshold.value;
	}

	//on equator - y is 0, near poles +-1 (convert it to +1) ->
	//we get percentual position in south and noth latitudes
	double yy = MyMath::MathUtils::MapRange<double>(
		0, rawData.GetHeight(),
		-1, 1,
		y);
	yy = std::abs(yy);

	return smallAreaThreshold.value * (latCorrectionFactor * yy + 1);
}

/// <summary>
/// Initialize mask used for extrema finding
/// Mask is square filled with symmetrical pattern
/// -1 -1 0 +1 +1
/// 
/// and
/// 
/// -1 
/// -1
///  0
/// +1
/// +1 
/// </summary>
/// <param name="maskX"></param>
/// <param name="maskY"></param>
void PressureFinder::InitExtremaMask(ConvolutionFilter2d & maskX,
	ConvolutionFilter2d & maskY)
{
	maskX = ConvolutionFilter2d::CreateSquare(maskRadius, 0);
	for (int y = 0; y < maskX.GetSize(); y++)
	{
		for (int x = 0; x < maskX.GetSize() / 2; x++)
		{
			maskX.GetValueAt(x, y)[0] = -1;
		}
	}
	for (int y = 0; y < maskX.GetSize(); y++)
	{
		for (int x = maskX.GetSize() / 2 + 1; x < maskX.GetSize(); x++)
		{
			maskX.GetValueAt(x, y)[0] = +1;
		}
	}

	//=========================================================================

	maskY = ConvolutionFilter2d::CreateSquare(maskRadius, 0);
	for (int y = 0; y < maskY.GetSize() / 2; y++)
	{
		for (int x = 0; x < maskY.GetSize(); x++)
		{
			maskY.GetValueAt(x, y)[0] = -1;
		}
	}
	for (int y = maskX.GetSize() / 2 + 1; y < maskY.GetSize(); y++)
	{
		for (int x = 0; x < maskY.GetSize(); x++)
		{
			maskY.GetValueAt(x, y)[0] = +1;
		}
	}	
}

/// <summary>
/// Calculate derivatives dx and dy of the input data
/// From derivatives, calculate signum values that are returned
/// </summary>
/// <returns></returns>
PressureFinder::Derivatives PressureFinder::CreateDerivatives()
{
	auto fDx = ConvolutionFilterSeparable::CreateSobelX();
	//auto fDx = ConvolutionFilter2d::CreateSobelX();
	auto dx = rawData.template CreateEmpty<float>();
	ImageFilters<float>::RunConvolutionFilter(rawData, dx, fDx);

	auto fDy = ConvolutionFilterSeparable::CreateSobelY();
	//auto fDy = ConvolutionFilter2d::CreateSobelY();
	auto dy = rawData.template CreateEmpty<float>();
	ImageFilters<float>::RunConvolutionFilter(rawData, dy, fDy);

	//=========================================================================


	NeighborhoodKernel k = NeighborhoodKernel::CreateFromRadius(maskRadius);

	Derivatives d;

	d.signDx = std::vector<float>(rawData.GetPixelsCount(), 0);
	d.signDy = std::vector<float>(rawData.GetPixelsCount(), 0);
	for (size_t i = 0; i < d.signDx.size(); i++)
	{
		d.signDx[i] = MyMath::MathUtils::Signum(dx.GetPixelStart(i)[0]);
		d.signDy[i] = MyMath::MathUtils::Signum(dy.GetPixelStart(i)[0]);
	}

	return d;
}

/// <summary>
/// Find extremas centers for Low and High pressure
/// Extremas are located inside contours that were marked as extrem candidate
/// </summary>
void PressureFinder::FindExtrema()
{	
#ifndef DISABLE_DEBUG_OUTPUT
	Image2d<uint8_t> test2;
	if (PressureFinder::debugDirName != "")
	{
		test2 = ColorSpace::ConvertGrayToRgb(rawData.CreateAsMapped<uint8_t>(0, 255));
	}
#endif

	ConvolutionFilter2d maskX;
	ConvolutionFilter2d maskY;
	this->InitExtremaMask(maskX, maskY);

	Derivatives d = this->CreateDerivatives();

	NeighborhoodKernel k = NeighborhoodKernel::CreateFromRadius(maskRadius);

	
	BVH<PressureExtrema *> bvh;	
	for (PressureExtrema & c : extremas)
	{				
		bvh.Add(&c, c.contour.c.GetBoundingBox());
	}	
	bvh.Build();
	

	for (auto type : { 
		PressureExtrema::PressureType::LOW, 
		PressureExtrema::PressureType::HIGH 
		})
	{
		MY_LOG_INFO("Processing: %s", (type == PressureExtrema::PressureType::LOW) ? "low" : "high");

		PressureExtrema * pxExtrema = nullptr;

		for (size_t i = 0; i < rawData.GetPixelsCount(); i++)
		{
			int xi, yi;
			rawData.GetPositionFromIndex(i, xi, yi);

			//check last extrema if the current pixel is inside
			//if yes - hooray, we dont have to test all contours
			if ((pxExtrema != nullptr) &&
				(pxExtrema->contour.c.IsPointInside({ xi, yi }) == false))
			{
				pxExtrema = nullptr;
			}


			if (pxExtrema == nullptr)
			{
				
				//point is not inside last extrema
				//=> test them all
				auto res = bvh.GetPointInside<Vector2d>({ xi, yi });
				for (auto r : res)
				{					
					if (r->contour.c.IsPointInside({ xi, yi }))
					{
						pxExtrema = r;
						break;
					}
				}
				
				/*
				//brute-force version of above code
				//test all extremas				
				for (PressureExtrema & c : extremas)
				{					
					if (c.contour.c.IsPointInside({ xi, yi }))
					{
						pxExtrema = &c;
						break;
					}
				}								
				*/
			}


			if (pxExtrema == nullptr)
			{
				//pixel is not inside extrema
				//do not process further
				continue;
			}

			int yy_start, yy_stop;
			int xx_start, xx_stop;
			ImageUtils::CalcKernelStartStop(xi, rawData.GetWidth() - 1, k, xx_start, xx_stop);
			ImageUtils::CalcKernelStartStop(yi, rawData.GetHeight() - 1, k, yy_start, yy_stop);


			int nonZeroArea[2] = { k.GetPixelArea() , k.GetPixelArea() };
			int isOk[2] = { 0, 0 };
			ConvolutionFilter2d * masks[2] = { &maskX, &maskY };
			std::vector<float> * imgD[2] = { &d.signDx, &d.signDy };

			bool extrem[2] = { false, false };

			//iterate pixel neighborhood given by mask
			//we run it 2x -> for dX and dY
			for (int j = 0; j < 2; j++)
			{										
				for (int yy = yy_start; yy <= yy_stop; yy++)
				{
					for (int xx = xx_start; xx <= xx_stop; xx++)
					{
						float mv = *(masks[j]->GetValueAt(xx - xx_start, yy - yy_start));
						if (mv != 0)
						{
							float tmp = (*imgD[j])[rawData.GetIndexFromPosition(xx, yy)];

							//test only pixels that are non zero in mask
							//in mask - zero pixels are at the center of the mask
							//Eg: for extrema min in X direction							
							//- - 0 + +

							//compare signum of derivative with mask signum
							if (tmp == mv)
							{
								//the signs are same - OK value
								isOk[j]++;
							}
							else if (tmp == 0.0f)
							{
								//if image derivation signum value is zero - decrease number of area pixels
								//we use number of area pixels to say how many pixels of this count
								//is OK value
								nonZeroArea[j]--;
							}
						}
					}
				}				
				
				//test if pixel [xi, yi] is extrema in direction X and Y
				//pixel is considered extrema if ratio of neighborhood pixels
				//that correspond to mask against "incorrect" pixels is above extremaCorrectRatio
				//(i.e. extremaCorrectRatio of area pixels is correct)
				extrem[j] = (((float)isOk[j] / nonZeroArea[j]) > extremaCorrectRatio);
				if (extrem[j] == false)
				{
					break;
				}
			}

			if (extrem[0] && extrem[1])
			{
				//there is correct xtrema in X and Y direction -> pixel is
				//function local extrema -> add pixel to the extrema
				pxExtrema->AddPixel(xi, yi);

#ifndef DISABLE_DEBUG_OUTPUT
				if (PressureFinder::debugDirName != "")
				{

					if (type == PressureExtrema::PressureType::HIGH)
					{
						test2[i][0] = 255;
						test2[i][1] = 0;
						test2[i][2] = 0;
					}
					else
					{
						test2[i][0] = 0;
						test2[i][1] = 0;
						test2[i][2] = 255;
					}
				}
#endif
			}
		}


		for (auto & c : extremas)
		{			
			c.CreateAreas(type, smallAreaThreshold);
		}
		
		//====== optional step ====
		//for each extrema center find the closest center in another
		//extrema

		if (minCentersDistanceThresholdPixel > 0)
		{
			float threshold2 = static_cast<float>(minCentersDistanceThresholdPixel * minCentersDistanceThresholdPixel);

			for (size_t i = 0; i < extremas.size(); i++)
			{
				auto& c1 = extremas[i];
				if (c1.IsValid() == false)
				{
					continue;
				}

				for (size_t j = i + 1; j < extremas.size(); j++)
				{
					auto& c2 = extremas[j];
					if (c2.IsValid() == false)
					{
						continue;
					}

					float dist = c1.center.DistanceSquared(c2.center);

					if ((dist != 0) && (dist < threshold2))
					{
						//extremas are too close
						//remove the smaller one
						if (c1.GetSizePixels() < c2.GetSizePixels())
						{
							c1.Clear();
						}
						else
						{
							c2.Clear();
						}
					}

				}
			}
		}
		//================================

		//we need to clear all pixels from extremas
		//because in the next iteration, HIGH pressure will be processed
		//the same area can have low / high pressure
		//and at the end, we choose the one with larger area
		for (auto & c : extremas)
		{
			c.ClearPixels();
		}


		maskX.Negate();
		maskY.Negate();
	}

#ifndef DISABLE_DEBUG_OUTPUT

	if (PressureFinder::debugDirName != "")
	{
		uint8_t col[3] = { 255, 255, 255 };

		for (auto& c : contours)
		{
			for (size_t i = 0; i < c.c.Size() - 1; i++)
			{
				auto p1 = c.c[i];
				auto p2 = c.c[i + 1];

#ifdef HAVE_OPENCV
				test2.RunOpenCV([&](auto img, auto m) {
					cv::line(m,
						cv::Point(int(p1.x), int(p1.y)),
						cv::Point(int(p2.x), int(p2.y)),
						cv::Scalar(col[0], col[1], col[2]),
						1);
					});
#else
				ImageUtils::DrawLine(test2, col,
					int(p1.x), int(p1.y),
					int(p2.x), int(p2.y));
#endif
			}
		}

		col[0] = 255;
		col[1] = 0;
		col[2] = 0;

		for (auto& c : extremas)
		{
			for (size_t i = 0; i < c.contour.c.Size() - 1; i++)
			{
				auto p1 = c.contour.c[i];
				auto p2 = c.contour.c[i + 1];


#ifdef HAVE_OPENCV
				test2.RunOpenCV([&](auto img, auto m) {
					cv::line(m,
						cv::Point(int(p1.x), int(p1.y)),
						cv::Point(int(p2.x), int(p2.y)),
						cv::Scalar(col[0], col[1], col[2]),
						2);
					});
#else
				ImageUtils::DrawLine(test2, col,
					int(p1.x), int(p1.y),
					int(p2.x), int(p2.y));
#endif
			}

		}

		std::string output = debugDirName;
		output += debugName;
		output += "_candidate_areas.png";

		test2.Save(output.c_str());
	}
#endif


}

void PressureFinder::RenderContours(const std::vector<Contour>& c) const
{
	if (PressureFinder::debugDirName == "")
	{
		return;
	}

	Image2d<uint8_t> pres = rawData.CreateAsMapped<uint8_t>(0, 255);
	Image2d<uint8_t> test = ColorSpace::ConvertGrayToRgb(pres);

	uint8_t col[3];

	for (auto& c : contours)
	{
		col[0] = (c.isExtrema) ? 255 : 255;
		col[1] = (c.isExtrema) ? 0 : 255;
		col[2] = (c.isExtrema) ? 0 : 255;
		int thick = (c.isExtrema) ? 2 : 1;



		for (size_t i = 0; i < c.c.Size() - 1; i++)
		{
			auto p1 = c.c[i];
			auto p2 = c.c[i + 1];

#ifdef HAVE_OPENCV
			test.RunOpenCV([&](auto img, auto m) {
				cv::line(m,
					cv::Point(int(p1.x), int(p1.y)),
					cv::Point(int(p2.x), int(p2.y)),
					cv::Scalar(col[0], col[1], col[2]),
					thick);
				});
#else
			ImageUtils::DrawLine(test, col,
				int(p1.x), int(p1.y),
				int(p2.x), int(p2.y));
#endif
		}
	}

	std::string output = debugDirName;
	output += debugName;
	output += "_all_contours_with_most_inner_highlighted.png";
	test.Save(output.c_str());
}

/// <summary>
/// Debug output if OpenCV is present during project build
/// </summary>
void PressureFinder::RenderContours()
{
	if (PressureFinder::debugDirName == "")
	{
		return;
	}

	Image2d<uint8_t> pres = rawData.CreateAsMapped<uint8_t>(0, 255);

	std::string output = PressureFinder::debugDirName;
	output += debugName;
	output += "_input_data.png";
	pres.Save(output.c_str());

	Image2d<uint8_t> test = ColorSpace::ConvertGrayToRgb(pres);

	uint8_t contourColor[3] = { 255, 255, 255 };
	uint8_t extremaColor[3] = { 255, 0, 0 };

	for (auto& c : contours)
	{
		for (size_t i = 0; i < c.c.Size() - 1; i++)
		{
			auto p1 = c.c[i];
			auto p2 = c.c[i + 1];

#ifdef HAVE_OPENCV
			test.RunOpenCV([&](auto img, auto m) {
				cv::line(m,
					cv::Point(int(p1.x), int(p1.y)),
					cv::Point(int(p2.x), int(p2.y)),
					cv::Scalar(255, 255, 255),
					1);
				});
#else
			ImageUtils::DrawLine(test, contourColor,
				int(p1.x), int(p1.y),
				int(p2.x), int(p2.y));
#endif

		}

		//lineIndex += 20;
	}

	for (auto& c : extremas)
	{
		for (size_t i = 0; i < c.contour.c.Size() - 1; i++)
		{
			auto p1 = c.contour.c[i];
			auto p2 = c.contour.c[i + 1];


#ifdef HAVE_OPENCV
			test.RunOpenCV([&](auto img, auto m) {
				cv::line(m,
					cv::Point(int(p1.x), int(p1.y)),
					cv::Point(int(p2.x), int(p2.y)),
					cv::Scalar(255, 0, 0),
					2);
				});
#else
			ImageUtils::DrawLine(test, extremaColor,
				int(p1.x), int(p1.y),
				int(p2.x), int(p2.y));
#endif
		}

	}

	output = debugDirName;
	output += debugName;
	output += "_extrema_candidate_contours_highligted.png";
	test.Save(output.c_str());


	test = ColorSpace::ConvertGrayToRgb(pres);


	for (auto& c : extremas)
	{
#ifdef HAVE_OPENCV
		test.RunOpenCV([&](auto img, auto mat) {

			cv::Scalar col = (c.type == PressureExtrema::PressureType::HIGH) ?
				cv::Scalar(255, 0, 0) :
				cv::Scalar(0, 0, 255);

			cv::circle(mat,
				cv::Point(int(c.center.x), int(c.center.y)),
				5, col, -5);

			});
#else
		uint8_t col[3] = { 0,0,0 };
		if (c.type == PressureExtrema::PressureType::HIGH)
		{
			col[0] = 255;
		}
		else
		{
			col[2] = 255;
		}

		ImageUtils::DrawLine(test, col,
			int(c.center.x - 5), int(c.center.y),
			int(c.center.x + 5), int(c.center.y));

		ImageUtils::DrawLine(test, col,
			int(c.center.x), int(c.center.y - 5),
			int(c.center.x), int(c.center.y + 5));
#endif		
	}



	output = debugDirName;
	output += debugName;
	output += "_extrema_location.png";
	test.Save(output.c_str());
}

/// <summary>
/// Get saver class instance
/// </summary>
/// <returns></returns>
PressureFinderSaver PressureFinder::GetSaver()
{
	return PressureFinderSaver(this);
}

//===============================================================================

PressureFinderSaver::PressureFinderSaver(const PressureFinder * pf) : 
	pf(pf)
{
}

PressureFinderSaver::~PressureFinderSaver()
{
}

void PressureFinderSaver::Save(const char * modelName, struct tm t,
	const char * dbFileName,
	std::function<void(int, int, Latitude &, Longitude &)> positionToGps)
{	
	

}

//=================================================================
//template specialization declaration

template double PressureFinder::CalcAreaProj(const Contour& c, Projections::Equirectangular* proj);
template double PressureFinder::CalcAreaProj(const Contour& c, Projections::Mercator* proj);
template double PressureFinder::CalcAreaProj(const Contour& c, Projections::LambertConic* proj);