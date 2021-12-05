#include "./PressureExtrema.h"

#include <stack>

#include "./MapProjections/Projections/Equirectangular.h"
#include "./MapProjections/Projections/Mercator.h"
#include "./MapProjections/Projections/LambertConic.h"

#include "./RasterData/Image2d.h"

#include "./Utils/Logger.h"

PressureExtrema::PressureExtrema(const Contour & c) :
	contour(c),
	bestPressure(std::nullopt),
	type(PressureType::NONE),
	center(MyMath::Vector2f::Infinity())
{
	this->ClearPixels();
}

/// <summary>
/// Test if extrema is valid (it has assigned best pressure area)
/// </summary>
/// <returns></returns>
bool PressureExtrema::IsValid() const
{
	return this->bestPressure.has_value();
}

/// <summary>
/// Add pixel to extrema
/// It will also update extremas AABB
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
void PressureExtrema::AddPixel(int x, int y)
{
	if (x < minX) minX = x;
	if (x > maxX) maxX = x;

	if (y < minY) minY = y;
	if (y > maxY) maxY = y;

	px.emplace_back(x, y);
}


Projections::Coordinate PressureExtrema::GetCenterGps(Projections::IProjectionInfo* proj) const
{
	if (auto eq = dynamic_cast<Projections::Equirectangular*>(proj))
	{
		return eq->ProjectInverse<float>({ center.x, center.y });
	}
	else if (auto me = dynamic_cast<Projections::Mercator*>(proj))
	{
		return me->ProjectInverse<float>({ center.x, center.y });
	}
	else if (auto lc = dynamic_cast<Projections::LambertConic*>(proj))
	{
		return lc->ProjectInverse<float>({ center.x, center.y });
	}

	MY_LOG_ERROR("Unknown projection");
	return Projections::Coordinate();
}

/// <summary>
/// Get center of AABB
/// </summary>
/// <param name="cx"></param>
/// <param name="cy"></param>
void PressureExtrema::GetAabbCenterIndex(int& cx, int& cy) const
{
	cx = minX + (maxX - minX) / 2;
	cy = minY + (maxY - minY) / 2;
}
int PressureExtrema::GetWidth() const
{
	return ((maxX - minX) + 1);
}

int PressureExtrema::GetHeight() const
{
	return ((maxY - minY) + 1);
}

/// <summary>
/// Get pixel area based on AABB
/// </summary>
/// <returns></returns>
double PressureExtrema::GetSizePixels() const
{
	return ((maxX - minX) + 1) * ((maxY - minY) + 1);
}

double PressureExtrema::GetSizeKm2(Projections::IProjectionInfo* proj) const
{
	if (auto eq = dynamic_cast<Projections::Equirectangular*>(proj))
	{
		std::vector<Projections::Pixel<float>> tmp;
		tmp.push_back({ float(minX), float(minY) });
		tmp.push_back({ float(maxX), float(minY) });
		tmp.push_back({ float(maxX), float(maxY) });
		tmp.push_back({ float(minX), float(maxY) });

		return Projections::ProjectionUtils::CalcArea(tmp, eq) / (1000.0 * 1000.0);
	}
	else if (auto me = dynamic_cast<Projections::Mercator*>(proj))
	{
		std::vector<Projections::Pixel<float>> tmp;
		tmp.push_back({ float(minX), float(minY) });
		tmp.push_back({ float(maxX), float(minY) });
		tmp.push_back({ float(maxX), float(maxY) });
		tmp.push_back({ float(minX), float(maxY) });

		return Projections::ProjectionUtils::CalcArea(tmp, me) / (1000.0 * 1000.0);
	}
	else if (auto lc = dynamic_cast<Projections::LambertConic*>(proj))
	{
		std::vector<Projections::Pixel<float>> tmp;
		for (size_t i = 0; i < px.size(); i++)
		{
			tmp.push_back({ float(px[i].x), float(px[i].y) });
		}

		return Projections::ProjectionUtils::CalcArea(tmp, lc) / (1000.0 * 1000.0);
	}

	MY_LOG_ERROR("Unknown projection - return area in pixels");
	return GetSizePixels();
}


/// <summary>
/// Clear enrire area
/// remove all pixels
/// </summary>
void PressureExtrema::ClearPixels()
{
	px.clear();

	minX = std::numeric_limits<int>::max();
	minY = std::numeric_limits<int>::max();
	maxX = 0;
	maxY = 0;
}

/// <summary>
/// Clear enrire area
/// remove all pixels
/// and remove best detected pressure
/// </summary>
void PressureExtrema::Clear()
{
	this->ClearPixels();

	this->type = PressureType::NONE;
	this->bestPressure = std::nullopt;
	this->center = MyMath::Vector2f::Infinity();
	this->centerValue = 0;
}


bool PressureExtrema::HasPixels() const
{
	return px.size() != 0;
}

bool PressureExtrema::CheckArea(const AreaThreshold& areaThreshold) const
{
	double area = 0;

	if (areaThreshold.unit == AreaThreshold::Unit::Km2)
	{
		area = this->GetSizeKm2(areaThreshold.proj);
	}
	else
	{
		area = this->GetSizePixels();
	}

	if (area <= areaThreshold.value)
	{
		//AABB area is too small (under threshold)		
		return false;
	}

	return true;
}

/// <summary>
/// Calculate center of extrema given by contour
/// We use pixels that were assigned to this contour 
/// The pixels may create multiple sub-areas
/// From the sub-areas, we use the largest one
/// </summary>
/// <param name="type"></param>
void PressureExtrema::InitCenter(PressureType type, const Image2d<float>& rawData)
{
	if (px.size() == 0)
	{
		//no pixels assigned to this extrema
		return;
	}

	//all pixels inside the contour can be "disconnected"
	//there can be multiple "components"
	//we use flood-fill to get connected components	
	//there is always at least one area
	std::list<PressureArea> areas = this->CreateAreasInner(rawData);

	for (auto& a : areas)
	{
		//take the larger area
		if ((bestPressure.has_value() == false) || (a.pixels.size() > bestPressure->pixels.size()))
		{
			bestPressure = a;

			this->type = type;
			if (type == PressureType::LOW)
			{
				this->center = MyMath::Vector2(bestPressure->minValuePx.x, bestPressure->minValuePx.y);
				this->centerValue = bestPressure->minValue;
			}
			else if (type == PressureType::HIGH)
			{
				this->center = MyMath::Vector2(bestPressure->maxValuePx.x, bestPressure->maxValuePx.y);
				this->centerValue = bestPressure->maxValue;
			}
			else
			{
				int cx, cy;
				bestPressure->GetAabbCenterIndex(cx, cy);
				this->center = MyMath::Vector2(cx, cy);
				this->centerValue = *rawData.GetPixelStart(cx, cy);
			}
		}
	}
}


void PressureExtrema::GetPositionFromIndex(int index, int & x, int & y) const
{
	int w = this->GetWidth();
	x = (index % w);
	y = (index / w);
}

int PressureExtrema::GetIndexFromPosition(int x, int y) const
{
	return x + y * this->GetWidth();
}

/// <summary>
/// Test if pixel at position [ind] is already visited or 
/// if it is pressure index in list of extrema pixels
/// </summary>
/// <param name="ind"></param>
/// <param name="pixels"></param>
/// <param name="result"></param>
/// <returns></returns>
bool PressureExtrema::TestAreaPixel(size_t ind, std::unordered_set<int> & pixels,
	std::unordered_map<int, PressureArea *> & result)
{
	if (result.find(ind) != result.end())
	{
		//already visited
		return false;
	}

	if (pixels.find(ind) == pixels.end())
	{
		//no pressure here
		return false;
	}

	return true;
}

/// <summary>
/// Flood-fill of all extrema pixels to create
/// connected components
/// Each component is single pressure area within the extrema
/// </summary>
/// <returns></returns>
std::list<PressureExtrema::PressureArea> PressureExtrema::CreateAreasInner(const Image2d<float>& rawData)
{
	//https://lodev.org/cgtutor/floodfill.html

	//printf("Building areas\n");

	//create list of pixel aligned to [0, 0] ... [AABB_w, AABB_h]
	std::unordered_set<int> pixels;
	for (const auto& p : px)
	{
		pixels.insert(this->GetIndexFromPosition(p.x - minX, p.y - minY));
	}

	std::list<PressureArea> fa; //list of single pressure areaa

	std::unordered_map<int, PressureArea*> result; //map of pixels to the area
													//[pixel index] = pointer to PressureArea to which pixel belongs

	int w = this->GetWidth();

	//iterate virtual image
	for (int y = minY; y <= maxY; y++)
	{
		for (int x = minX; x <= maxX; x++)
		{
			int i = (x - minX) + (y - minY) * w;

			if (TestAreaPixel(i, pixels, result) == false)
			{
				continue;
			}

			ImageUtils::Pixel px;
			this->GetPositionFromIndex(i, px.x, px.y);

			//create bounding box in real pixels 
			//(translate pixels back to their original position)
			PressureArea* f = &(fa.emplace_back());
			f->minX = (px.x + minX);
			f->minY = (px.y + minY);
			f->maxX = (px.x + minX);
			f->maxY = (px.y + minY);

			f->pixels.reserve(100);

			//flood the pixels neighborhood

			int x1;
			bool spanAbove, spanBelow;

			std::stack<ImageUtils::Pixel> q;
			q.push(px);

			while (q.empty() == false)
			{
				ImageUtils::Pixel px = q.top();
				q.pop();

				x1 = px.x;
				int ind = px.y * w + x1;


				//test if a pixel is part of the currently flooded area		
				while ((x1 >= 0) && (TestAreaPixel(ind, pixels, result)))
				{
					x1--;
					ind--;
				}

				x1++;
				ind++;

				spanAbove = spanBelow = false;

				//test if a pixel is part of the currently flooded area
				while ((x1 < w) && (TestAreaPixel(ind, pixels, result)))
				{
					{
						//recreate pixel [xx, yy] in the original image
						int xx, yy;
						this->GetPositionFromIndex(ind, xx, yy);

						xx += minX;
						yy += minY;

						if (xx < f->minX) f->minX = xx;
						else if (xx > f->maxX) f->maxX = xx;

						if (yy < f->minY) f->minY = yy;
						else if (yy > f->maxY) f->maxY = yy;

						//add pixel to area				
						f->pixels.emplace_back(xx, yy);

						float val = *rawData.GetPixelStart(xx, yy);
						if (val > f->maxValue)
						{
							f->maxValue = val;
							f->maxValuePx.x = xx;
							f->maxValuePx.y = yy;
						}
						if (val < f->minValue)
						{
							f->minValue = val;
							f->minValuePx.x = xx;
							f->minValuePx.y = yy;
						}

						result[ind] = f;
					}

					if (px.y > 0)
					{
						//test if a pixel is part of the currently flooded area
						bool aboveRes = TestAreaPixel(ind - w, pixels, result);

						if ((!spanAbove) && (aboveRes))
						{
							q.push({ x1, px.y - 1 });
							spanAbove = true;
						}
						else if ((spanAbove) && (!aboveRes))
						{
							spanAbove = false;
						}
					}

					if (px.y < this->GetHeight() - 1)
					{
						//test if a pixel is part of the currently flooded area
						bool belowRes = TestAreaPixel(ind + w, pixels, result);

						if ((!spanBelow) && (belowRes))
						{
							q.push({ x1, px.y + 1 });
							spanBelow = true;
						}
						else if ((spanBelow) && (!belowRes))
						{
							spanBelow = false;
						}
					}

					x1++;
					ind++;
				}
			}
		}
	}


	return fa;
}

