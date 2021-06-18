#include "./ImageUtils.h"

#include <stdlib.h>
#include <algorithm>
#include <stack>
#include <queue>
#include <thread>

#ifdef ENABLE_SIMD
#	include <immintrin.h>
#endif

#include "../Utils/Logger.h"

#include "./Image2d.h"
#include "./ImageFilters.h"
#include "./ImageGausFilter.h"

//=================================================================================================
// Neighbothood class
//==============================================================================================

/// <summary>
/// Create kernel from overal size
/// Overal size is total length of pixels
/// M x M
/// Radius is M / 2
/// 
/// Note: most of the time, this value should be odd (3, 5, 7 ...) 
/// </summary>
/// <param name="size"></param>
/// <returns></returns>
NeighborhoodKernel NeighborhoodKernel::CreateFromSize(int size)
{
	return NeighborhoodKernel(size / 2);
}

/// <summary>
/// Create kernel from radius
/// It is number of pixels in a ring around the center
/// Size of kernel is r * 2 + 1
/// </summary>
/// <param name="r"></param>
/// <returns></returns>
NeighborhoodKernel NeighborhoodKernel::CreateFromRadius(int r)
{
	return NeighborhoodKernel(r);
}

/// <summary>
/// private ctor
/// instance are created via Factory methods
/// </summary>
/// <param name="radius"></param>
NeighborhoodKernel::NeighborhoodKernel(int radius) :
	radius(radius)
{
}

/// <summary>
/// Get size of filter
/// </summary>
/// <returns></returns>
int NeighborhoodKernel::GetSize() const
{
	return radius * 2 + 1;
}

/// <summary>
/// Get number of pixels the kernel covers
/// </summary>
/// <returns></returns>
int NeighborhoodKernel::GetPixelArea() const
{
	int s = this->GetSize();
	return s * s;
}

//=================================================================================================
// Drawing
//==============================================================================================

/// <summary>
/// Create random RGB color
/// Use rand(). It must be inited before
/// </summary>
/// <returns></returns>
std::array<uint8_t, 3> ImageUtils::CreateRandomColor()
{	
	return { uint8_t(rand() & 255), uint8_t(rand() & 255), uint8_t(rand() & 255) };
}

/// <summary>
/// Codes for Cohen-Sutherland
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
int ImageUtils::ComputeOutCode(int x, int y, int w, int h)
{
	int code = INSIDE;

	if (x < 0) code |= LEFT;
	else if (x >= w) code |= RIGHT;

	if (y < 0) code |= BOTTOM;
	else if (y >= h) code |= TOP;

	return code;
}


/// <summary>
/// Draw simple line (no anti-aliasing) to image with simple Bresenham algorithm
/// Line if from: [x0, y0] -> [x1, y1] 
/// If line is outside image bounds, it is clamped with Cohen-Sutherland
/// </summary>
/// <param name="input"></param>
/// <param name="value"></param>
/// <param name="x0"></param>
/// <param name="y0"></param>
/// <param name="x1"></param>
/// <param name="y1"></param>
template <typename T>
void ImageUtils::DrawLine(Image2d<T> & input, const T * value,
	int x0, int y0, int x1, int y1)
{
	// compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
	int outcode0 = ImageUtils::ComputeOutCode(x0, y0, input.GetWidth(), input.GetHeight());
	int outcode1 = ImageUtils::ComputeOutCode(x1, y1, input.GetWidth(), input.GetHeight());
	bool accept = false;

	double xmin = 0;
	double xmax = input.GetWidth() - 1;

	double ymin = 0;
	double ymax = input.GetHeight() - 1;

	while (true)
	{
		if (!(outcode0 | outcode1))
		{
			// Bitwise OR is 0. Trivially accept and get out of loop
			accept = true;
			break;
		}
		else if (outcode0 & outcode1)
		{
			// Bitwise AND is not 0. (implies both end points are in the same region outside the window). Reject and get out of loop
			break;
		}
		else
		{
			// failed both tests, so calculate the line segment to clip
			// from an outside point to an intersection with clip edge
			double x = 0;
			double y = 0;

			// At least one endpoint is outside the clip rectangle; pick it.
			int outcodeOut = outcode0 ? outcode0 : outcode1;

			// Now find the intersection point;
			// use formulas y = y0 + slope * (x - x0), x = x0 + (1 / slope) * (y - y0)
			if (outcodeOut & TOP) {           // point is above the clip rectangle
				x = x0 + (x1 - x0) * (ymax - y0) / (y1 - y0);
				y = ymax;
			}
			else if (outcodeOut & BOTTOM) { // point is below the clip rectangle
				x = x0 + (x1 - x0) * (ymin - y0) / (y1 - y0);
				y = ymin;
			}
			else if (outcodeOut & RIGHT) {  // point is to the right of clip rectangle
				y = y0 + (y1 - y0) * (xmax - x0) / (x1 - x0);
				x = xmax;
			}
			else if (outcodeOut & LEFT) {   // point is to the left of clip rectangle
				y = y0 + (y1 - y0) * (xmin - x0) / (x1 - x0);
				x = xmin;
			}

			// Now we move outside point to intersection point to clip
			// and get ready for next pass.
			if (outcodeOut == outcode0)
			{
				x0 = static_cast<int>(x);
				y0 = static_cast<int>(y);
				outcode0 = ImageUtils::ComputeOutCode(x0, y0, input.GetWidth(), input.GetHeight());
			}
			else
			{
				x1 = static_cast<int>(x);
				y1 = static_cast<int>(y);
				outcode1 = ImageUtils::ComputeOutCode(x1, y1, input.GetWidth(), input.GetHeight());
			}
		}
	}

	if (accept == false)
	{
		return;
	}

	ImageUtils::ProcessLinePixels(x0, y0, x1, y1,
		[&](int x, int y) {

		auto * val = input.GetPixelStart(x, y);
		for (size_t c = 0; c < input.GetChannelsCount(); c++)
		{
			val[c] = value[c];
		}

	});
}

/// <summary>
/// Iterate pixels on simple line (no anti-aliasing) with simple Bresenham algorithm
/// Line if from: [x0, y0] -> [x1, y1] 
/// 
/// For each line pixel, pixelCallback(x, y) is called
/// </summary>
/// <param name="x0"></param>
/// <param name="y0"></param>
/// <param name="x1"></param>
/// <param name="y1"></param>
/// <param name="pixelCallback"></param>
void ImageUtils::ProcessLinePixels(int x0, int y0, int x1, int y1,
	std::function<void(int x, int y)> pixelCallback)
{
			   
	int dx = std::abs(x1 - x0);
	int dy = std::abs(y1 - y0);
	int sx, sy, e2;


	(x0 < x1) ? sx = 1 : sx = -1;
	(y0 < y1) ? sy = 1 : sy = -1;
	int err = dx - dy;

	while (1)
	{		
		pixelCallback(x0, y0);
		
		if ((x0 == x1) && (y0 == y1))
		{
			break;
		}
		e2 = 2 * err;
		if (e2 > -dy)
		{
			err = err - dy;
			x0 = x0 + sx;
		}
		if (e2 < dx)
		{
			err = err + dx;
			y0 = y0 + sy;
		}
	}
}

/// <summary>
/// Queue-based flood fill with a valDifThreshold
/// value is filled if pixel difference against newVal is 
/// under valDifThreshold
/// </summary>
/// <param name="input"></param>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="newVal">new filling value</param>
/// <param name="valDifThreshold">difference between orig and current pixel value </param>
/// <param name="con">type of neighborhood connectivity</param>
template <typename T>
size_t ImageUtils::FloodFill(Image2d<T> & input, int x, int y,
	T newVal, T valDifThreshold,
	NeighborhoodKernel::Connectivity con)
{
	return FloodFill<T>(input, x, y, newVal,
		[&](T orig, T px, size_t ind) -> bool {
		return std::abs(px - orig) < valDifThreshold;
	}, con);
		
}

/// <summary>
/// Queue-based flood fill with a valDifThreshold
/// value is filled if pixel difference against newVal is 
/// under valDifThreshold
/// To test if value belongs to the area, 
/// sameAreaCallback(originalValue, currentValue, currenrValueIndex)
/// is called
/// </summary>
/// <param name="input"></param>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="newVal">new filling value</param>
/// <param name="sameAreaCallback">test callback if pixel is in same area</param>
/// <param name="con">type of neighborhood connectivity</param>
template <typename T>
size_t ImageUtils::FloodFill(Image2d<T> & input, int x, int y,
	T newVal, 
	std::function<bool(T, T, size_t)> sameAreaCallback,
	NeighborhoodKernel::Connectivity con)
{
	size_t filled = 0;

	Image2d<uint8_t> visited = input.template CreateEmpty<uint8_t>();

	size_t origInd = input.GetIndexFromPosition(x, y);
		
	T origVal = input.GetPixelStart(origInd)[0];
	
	std::queue<size_t> q;
	q.push(origInd);
	
	while (q.empty() == false)
	{
		size_t ind = q.front();
		q.pop();

		T * px = input.GetPixelStart(ind);
		uint8_t * vis = visited.GetPixelStart(ind);

		if ((vis[0] == 0) && (sameAreaCallback(origVal, px[0], ind)))
		{		
			vis[0] = 255;
			px[0] = newVal;
			filled++;

			input.GetPositionFromIndex(ind, x, y);


			int x1m = (x - 1 < 0) ? 0 : x - 1;
			int x1p = (x + 1 >= input.GetWidth()) ? input.GetWidth() - 1 : x + 1;

			int y1m = (y - 1 < 0) ? 0 : y - 1;
			int y1p = (y + 1 >= input.GetHeight()) ? input.GetHeight() - 1 : y + 1;

			if (con == NeighborhoodKernel::Connectivity::CON_8)
			{
				q.push(input.GetIndexFromPosition(x1m, y1m));
				q.push(input.GetIndexFromPosition(x1m, y));
				q.push(input.GetIndexFromPosition(x1m, y1p));

				q.push(input.GetIndexFromPosition(x, y1m));
				q.push(input.GetIndexFromPosition(x, y1p));

				q.push(input.GetIndexFromPosition(x1p, y1m));
				q.push(input.GetIndexFromPosition(x1p, y));
				q.push(input.GetIndexFromPosition(x1p, y1p));
			}
			else
			{
				q.push(input.GetIndexFromPosition(x1m, y));
				q.push(input.GetIndexFromPosition(x1p, y));
				q.push(input.GetIndexFromPosition(x, y1m));
				q.push(input.GetIndexFromPosition(x, y1p));
			}
		}
	}

	return filled;
}

//==============================================================================================
// Vectorization
//==============================================================================================

/// <summary>
/// Vectorize binary image using simple DFS-based technique
/// 
/// Lines shorter than minLength2 are removed
/// If joining lines is enabled, neighboring lines are join together if they share corner at 3x3 neighborhood
/// </summary>
/// <param name="input"></param>
/// <param name="bgValue">background value - not to be vectorized</param>
/// <param name="minLength2">minimal squared length of vectorized line (default: 0)</param>
/// <param name="joinLines">(default: false)</param>
/// <returns></returns>
template <typename T>
std::list<std::vector<ImageUtils::Pixel>> ImageUtils::Vectorize(const Image2d<T> & input, T bgValue, 
	double minLength2, bool joinLines)
{
	//create copy of input
	//in this, visited pixels are erased
	Image2d<T> inputCopy = input.CreateDeepCopy();

	std::list<std::vector<ImageUtils::Pixel>> lines;
	NeighborhoodKernel n3 = NeighborhoodKernel::CreateFromSize(3);

	for (size_t i = 0; i < inputCopy.GetPixelsCount(); i++)
	{
		if (inputCopy[i][0] == bgValue)
		{
			continue;
		}

		std::stack<size_t> st;
		inputCopy[i][0] = bgValue;
		st.push(i);

		std::vector<ImageUtils::Pixel> p;
		int sinceLastSplitCount = 0;
		double lineLen = 0;

		while (st.empty() == false)
		{
			size_t ii = st.top();
			st.pop();

			int addedCount = 0;
			int x, y;
			inputCopy.GetPositionFromIndex(ii, x, y);

			inputCopy.ForEachPixelInNeighborhood(x, y, n3,
				[&](T * val, int xx, int yy) {
				if (val[0] == bgValue)
				{
					return;
				}

				size_t neighIndex = inputCopy.GetIndexFromPosition(xx, yy);
				st.push(neighIndex);
				inputCopy[neighIndex][0] = bgValue;
				addedCount++;
			});

			if (p.size() > 1)
			{
				lineLen += (p.back().x - x) * (p.back().x - x) + (p.back().y - y) * (p.back().y - y);
			}

			p.emplace_back(x, y);
						
			if (addedCount == 0)
			{		
				ImageUtils::AddVectorizedLine(p, sinceLastSplitCount,
					lineLen, minLength2,
					lines);

				p.clear();
				lineLen = 0;
				sinceLastSplitCount = 0;
			}
			else if (addedCount == 1)
			{
				sinceLastSplitCount++;
			}
			else 
			{
				sinceLastSplitCount = 0;
			}

		}

		//last segment
		ImageUtils::AddVectorizedLine(p, sinceLastSplitCount,
			lineLen, minLength2,
			lines);		
	}
	
	if (joinLines)
	{
		ImageUtils::JoinVectorizedLines(input.GetWidth(), input.GetHeight(), lines);
	}
			
	return lines;
}

/// <summary>
/// Add vectorized part to the "global" list of all vectorized lines
/// </summary>
/// <param name="newLine"></param>
/// <param name="sinceLastSplitCount"></param>
/// <param name="lineLen"></param>
/// <param name="minLength2"></param>
/// <param name="dest"></param>
void ImageUtils::AddVectorizedLine(std::vector<ImageUtils::Pixel> & newLine,
	int sinceLastSplitCount,
	double lineLen,
	double minLength2,
	std::list<std::vector<Pixel>> & dest)
{
	//check if the first point of a new line corresponds to 
	//the 3x3 neighbohood of the last point of the last line
	if ((dest.size() > 0) && (newLine.size() > 0))
	{
		
		auto & ll = dest.back();
		int dx = std::abs(ll.back().x - newLine.front().x);
		int dy = std::abs(ll.back().y - newLine.front().y);
		if ((dx <= 1) && (dy <= 1))
		{
			//append line
			for (size_t r = 0; r < newLine.size(); r++)
			{
				ll.emplace_back(newLine[r].x, newLine[r].y);
			}

			newLine.clear();
			sinceLastSplitCount = std::numeric_limits<int>::max();
		}
		
	}

	if (sinceLastSplitCount < 5)
	{
		//line length since last "intersection" (split) is too short
		//probably a "dangling" branch - remove it
		for (int r = 0; r <= sinceLastSplitCount; r++)
		{
			newLine.pop_back();
		}
	}

	if ((newLine.size() > 1) && (lineLen >= minLength2))
	{
		dest.emplace_back(std::move(newLine));
	}
}

/// <summary>
/// Join vectorized lines
/// joining occurs if lines share corner or edge at 3x3 neighborhood
/// Empty lines are removed from lines list
/// </summary>
/// <param name="w"></param>
/// <param name="h"></param>
/// <param name="lines"></param>
void ImageUtils::JoinVectorizedLines(int w, int h, std::list<std::vector<Pixel>> & lines)
{
	//connect 

	NeighborhoodKernel k = NeighborhoodKernel::CreateFromRadius(1);
	int yy_start0, yy_stop0;
	int xx_start0, xx_stop0;

	int yy_start1, yy_stop1;
	int xx_start1, xx_stop1;

	std::vector<std::list<std::vector<Pixel>>::iterator> removed;

	std::list<std::vector<Pixel>>::iterator it = lines.begin();
	while (it != lines.end())
	{				
		if (it->size() < 2)
		{
			it++;
			continue;
		}

		ImageUtils::CalcKernelStartStop((*it).back().x, w - 1, k, xx_start1, xx_stop1);
		ImageUtils::CalcKernelStartStop((*it).back().y, h - 1, k, yy_start1, yy_stop1);


		ImageUtils::CalcKernelStartStop((*it)[0].x, w - 1, k, xx_start0, xx_stop0);
		ImageUtils::CalcKernelStartStop((*it)[0].y, h - 1, k, yy_start0, yy_stop0);

		std::list<std::vector<Pixel>>::iterator jt;
		for (jt = lines.begin(); jt != lines.end(); jt++)
		{
			if (it == jt) continue;
			if ((*jt).size() < 2) continue;

			int startX = (*jt)[0].x;
			int startY = (*jt)[0].y;

			int endX = (*jt).back().x;
			int endY = (*jt).back().y;
			
			
			for (int yy = yy_start0; yy <= yy_stop0; yy++)
			{
				for (int xx = xx_start0; xx <= xx_stop0; xx++)
				{
					if ((xx == endX) && (yy == endY))					
					{
						//join lines (start -> end)
						(*jt).insert((*jt).end(), (*it).begin(), (*it).end());
						(*it).clear();
						removed.push_back(it);
							
						it = lines.begin();
						goto joinend;
					}

					if ((xx == startX) && (yy == startY))
					{
						//join lines (start -> start .. reverse first line)
						std::reverse((*it).begin(), (*it).end());
						(*it).insert((*it).end(), (*jt).begin(), (*jt).end());
						(*jt).clear();
						removed.push_back(jt);

						it = lines.begin();
						goto joinend;
					}
				}
			}
			

			for (int yy = yy_start1; yy <= yy_stop1; yy++)
			{
				for (int xx = xx_start1; xx <= xx_stop1; xx++)
				{
					if ((xx == startX) && (yy == startY))
					{
						//join lines (end -> start)
						(*it).insert((*it).end(), (*jt).begin(), (*jt).end());
						(*jt).clear();
						removed.push_back(jt);

						it = lines.begin();
						goto joinend;
					}
					if ((xx == endX) && (yy == endY))
					{
						//join lines (end -> end .. reverse second line)
						std::reverse((*jt).begin(), (*jt).end());
						(*it).insert((*it).end(), (*jt).begin(), (*jt).end());
						(*jt).clear();
						removed.push_back(jt);

						it = lines.begin();
						goto joinend;
					}
				}
			}									
		}

		it++;
	joinend: ;
		
	}

	for (auto r : removed)
	{
		lines.erase(r);
	}
}

//==============================================================================================
// SDF
// Source: https://github.com/prideout/heman/blob/master/src/distance.c
//==============================================================================================

const float ImageUtils::INF = 1E20f;

template <typename T>
Image2d<float> ImageUtils::SDF(const Image2d<T> & input, T bgValue)
{
	Image2d<float> positive = input.template CreateEmpty<float>();
	Image2d<float> negative = input.template CreateEmpty<float>();
	
	const size_t size = input.GetPixelsCount();

	float * pptr = positive.GetData().data();
	float * nptr = negative.GetData().data();
	const T * sptr = input.GetData().data();
	
	for (size_t i = 0; i < size; ++i, ++sptr)
	{
		bool isBg = (*sptr == bgValue);
		*pptr++ = (isBg) ? INF : 0;
		*nptr++ = (isBg) ? 0 : INF;
	}

	std::thread t1([&] {
		SDFTransformToDistance(positive);
	});

	std::thread t2([&] {
		SDFTransformToDistance(negative);
	});
	t1.join();
	t2.join();

	const float inv = 1.0f / input.GetWidth();
	pptr = positive.GetData().data();
	nptr = negative.GetData().data();

	size_t size8 = 0;
#ifdef ENABLE_SIMD

	__m256 inv8 = _mm256_set1_ps(inv);

	size8 = size - size % 8;
	for (size_t i = 0; i < size8; i += 8)
	{
		__m256 a = _mm256_loadu_ps(pptr);
		__m256 b = _mm256_loadu_ps(nptr);
			
		a = _mm256_sqrt_ps(a);
		b = _mm256_sqrt_ps(b);

		a = _mm256_sub_ps(a, b);
		a = _mm256_mul_ps(a, inv8);

		_mm256_storeu_ps(pptr, a);

		pptr += 8;
		nptr += 8;
	}


#endif
	for (size_t i = size8; i < size; ++i, ++pptr, ++nptr) 
	{
		*pptr = (sqrt(*pptr) - sqrt(*nptr)) * inv;
	}

	return positive;
}


void ImageUtils::SDFTransformToDistance(Image2d<float> & img)
{		
	float* fh = (float *)malloc(img.GetHeight() * sizeof(float));
	
	int hwMax = std::max(img.GetHeight(), img.GetWidth());

	//allocate helper arrays to maximum of width or height
	//we allocate them once and use them for both pass
	float* zhw = (float *)malloc((hwMax + 1) * sizeof(float));	
	float* dhw = (float *)malloc(hwMax * sizeof(float));	
	
	uint16_t* hw = (uint16_t *)calloc(hwMax, sizeof(uint16_t)); //need to be zeroed
	

	/*
	float* zh = (float *)calloc((img.GetHeight() + 1), sizeof(float));
	float* zw = (float *)calloc((img.GetWidth() + 1), sizeof(float));
	
	float* dh = (float *)calloc(img.GetHeight(), sizeof(float));
	float* dw = (float *)calloc(img.GetWidth(), sizeof(float));

	uint16_t* hh = (uint16_t *)calloc(img.GetHeight(), sizeof(uint16_t));
	uint16_t* ww = (uint16_t *)calloc(img.GetWidth(), sizeof(uint16_t));	
	*/

	for (int x = 0; x < img.GetWidth(); ++x) 
	{		
		for (int y = 0; y < img.GetHeight(); ++y) 
		{				
			fh[y] = img.GetPixelStart(x, y)[0];
		}
		SDFCalcEdt(fh, dhw, zhw, hw, img.GetHeight());
		for (int y = 0; y < img.GetHeight(); ++y) 
		{			
			img.GetPixelStart(x, y)[0] = dhw[y];
		}
	}

	//need to be zeroed
	memset(hw, 0, hwMax * sizeof(uint16_t));

	for (int y = 0; y < img.GetHeight(); ++y) 
	{
		float * f = img.GetData().data() + (y * img.GetWidth());

		SDFCalcEdt(f, dhw, zhw, hw, img.GetWidth());
		for (int x = 0; x < img.GetWidth(); ++x) 
		{		
			img.GetPixelStart(x, y)[0] = dhw[x];
		}
	}


	free(fh);
	
	free(zhw);
	free(dhw);
	free(hw);

	/*
	free(dw);
	free(dh);
	free(zh);
	free(zw);
	free(ww);
	free(hh);
	*/
}

void ImageUtils::SDFCalcEdt(const float * f, float* d, float* z, uint16_t* w, int n)
{
	int k = 0;
	float s;
	w[0] = 0;
	z[0] = -INF;
	z[1] = +INF;

	for (int q = 1; q < n; ++q) 
	{
		uint16_t wk = w[k];
		float fq2 = f[q] + (q * q);
		s = (fq2 - (f[wk] + (wk * wk))) / (2 * (q - wk));
		
		//this usually runs 0 or 1 times
		while (s <= z[k]) 
		{					
			--k;
			wk = w[k];			
			s = (fq2 - (f[wk] + (wk * wk))) / (2 * (q - wk));
		}

		w[++k] = q;
		z[k] = s;
		z[k + 1] = +INF;
	}

	k = 0;
	for (int q = 0; q < n; ++q) 
	{
		while (z[k + 1] < q) 
		{
			++k;
		}
		d[q] = ((q - w[k]) * (q - w[k])) + f[w[k]];
	}
}

//==============================================================================================
// Utils
//==============================================================================================

/// <summary>
/// Helper method to calculate kernel start and end position in image
/// </summary>
/// <param name="v0"></param>
/// <param name="maxVal"></param>
/// <param name="kernelSize"></param>
/// <param name="start"></param>
/// <param name="end"></param>
void ImageUtils::CalcKernelStartStop(double v0, int maxVal,
	const NeighborhoodKernel & k,
	int & start, int & end)
{
	int c = static_cast<int>(std::floor(v0));
	return ImageUtils::CalcKernelStartStop(c, maxVal, k, start, end);
}

/// <summary>
/// Helper method to calculate kernel start and end position in image
/// </summary>
/// <param name="v0"></param>
/// <param name="maxVal"></param>
/// <param name="kernelSize"></param>
/// <param name="start"></param>
/// <param name="end"></param>
void ImageUtils::CalcKernelStartStop(int v0, int maxVal,
	const NeighborhoodKernel & k,
	int & start, int & end)
{		
	start = std::max(v0 - k.radius, 0);
	end = std::min(v0 + k.radius, maxVal);
}

/// <summary>
/// Get the heigts maps first and second derivative using Evans-Young method.
/// </summary>
/// <param name="input"></param>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="dx"></param>
/// <param name="dy"></param>
/// <param name="dxx"></param>
/// <param name="dyy"></param>
/// <param name="dxy"></param>
template <typename T>
void ImageUtils::CalcDerivatives(const Image2d<T> & input, int x, int y,
	float * dx, float * dy, float * dxx, float * dyy, float * dxy)
{
	const int x1m = (x < 1) ? (0) : (x - 1);
	const int x1p = (x + 1 >= input.GetWidth()) ? (input.GetWidth() - 1) : (x + 1);

	const int y1m = (y < 1) ? (0) : (y - 1);
	const int y1p = (y + 1 >= input.GetHeight()) ? (input.GetHeight() - 1) : (y + 1);
	
	const float z1 = (float)(input.GetPixelStart(x1m, y1m)[0]);
	const float z2 = (float)(input.GetPixelStart(x  , y1m)[0]);
	const float z3 = (float)(input.GetPixelStart(x1p, y1m)[0]);

	const float z4 = (float)(input.GetPixelStart(x1m, y)[0]);
	const float z5 = (float)(input.GetPixelStart(x  , y)[0]);
	const float z6 = (float)(input.GetPixelStart(x1p, y)[0]);

	const float z7 = (float)(input.GetPixelStart(x1m, y1p)[0]);
	const float z8 = (float)(input.GetPixelStart(x  , y1p)[0]);
	const float z9 = (float)(input.GetPixelStart(x1p, y1p)[0]);

	//p, q
	if (dx) *dx = (z3 + z6 + z9 - z1 - z4 - z7) / (6.0f);
	if (dy) *dy = (z1 + z2 + z3 - z7 - z8 - z9) / (6.0f);

	//r, t, s
	if (dxx) *dxx = (z1 + z3 + z4 + z6 + z7 + z9 - 2.0f * (z2 + z5 + z8)) / (3.0f);
	if (dyy) *dyy = (z1 + z2 + z3 + z7 + z8 + z9 - 2.0f * (z4 + z5 + z6)) / (3.0f);
	if (dxy) *dxy = (z3 + z7 - z1 - z9) / (4.0f);
}


//==============================================================================================
// Corner detection
//==============================================================================================

/// <summary>
/// Detect corners in image
/// return image with corner response function
/// bestCorners is optional output, that will fill array
/// with best corners
/// </summary>
/// <param name="input"></param>
/// <param name="info"></param>
/// <param name="bestCorners"></param>
/// <returns></returns>
template <typename T>
Image2d<float> ImageUtils::DetectCornersHarris(const Image2d<T> & input, const CornerDetection & info,
	std::vector<size_t> * bestCorners)
{
	Image2d<float> res = input.template CreateEmpty<float>();

	if(input.GetChannelsCount() != 1)
	{
		MY_LOG_ERROR("Corner Harris detector input must have only 1 channel");
		return res;
	}
		
	Image2d<float> dx, dy;

	if (input.GetPixelsCount() > 1048576)
	{
		//too many pixels
		//we will run dx and dy calculation in parallel
		//creation of threads is faster than 2 serial runs

		std::thread t1([&] {
			ConvolutionFilterSeparable fDx = ConvolutionFilterSeparable::CreateSobelX();
			dx = input.template CreateEmpty<float>();
			ImageFilters<T>::RunConvolutionFilter(input, dx, fDx);
		});
		std::thread t2([&] {
			ConvolutionFilterSeparable fDy = ConvolutionFilterSeparable::CreateSobelY();
			dy = input.template CreateEmpty<float>();
			ImageFilters<T>::RunConvolutionFilter(input, dy, fDy);
		});
		t1.join();
		t2.join();
	}
	else
	{
		ConvolutionFilterSeparable fDx = ConvolutionFilterSeparable::CreateSobelX();
		dx = input.template CreateEmpty<float>();
		ImageFilters<T>::RunConvolutionFilter(input, dx, fDx);

		ConvolutionFilterSeparable fDy = ConvolutionFilterSeparable::CreateSobelY();
		dy = input.template CreateEmpty<float>();
		ImageFilters<T>::RunConvolutionFilter(input, dy, fDy);
	}



	//Compute elements of the local structure matrix M
	auto mulCallback = [&](float * a, const float * b, size_t c) {
		a[0] *= b[0];		
	};
	
	
	Image2d<float> dxy = dx.CreateDeepCopy();
	dxy.Combine(dy, mulCallback); //dx * dy	
	dx.Combine(dx, mulCallback); //dx * dx
	dy.Combine(dy, mulCallback); //dy * dy

	Image2d<float> dx2 = std::move(dx);
	Image2d<float> dy2 = std::move(dy);


	//Blur each component of the structure matrix M
	if (input.GetPixelsCount() > 1048576)
	{
		std::thread t1([&] {
			dx2.RunGauss(info.gausRadius, info.gausWeight, info.gausStrength);
		});
		std::thread t2([&] {
			dy2.RunGauss(info.gausRadius, info.gausWeight, info.gausStrength);
		});
		std::thread t3([&] {
			dxy.RunGauss(info.gausRadius, info.gausWeight, info.gausStrength);
		});
		t1.join();
		t2.join();
		t3.join();
	}
	else
	{
		dx2.RunGauss(info.gausRadius, info.gausWeight, info.gausStrength);
		dy2.RunGauss(info.gausRadius, info.gausWeight, info.gausStrength);
		dxy.RunGauss(info.gausRadius, info.gausWeight, info.gausStrength);
	}
	
	//Compute the corner response function
	
	size_t len = input.GetPixelsCount();
	for (size_t index = 0; index < len; index++)
	{
		//calculate response function from M
		float a = dx2.GetData()[index];
		float b = dy2.GetData()[index];
		float c = dxy.GetData()[index];

		float q = info.cornerFunctionCallback(a, b, c, info.kappa);
		
		res.GetData()[index] = q;
	}
	
	
	if (bestCorners == nullptr)
	{
		return res;
	}

	
	if ((info.k.radius != 0) && (info.percentage != 0))
	{
		*bestCorners = ImageUtils::FindMaxValues(res, info.k, info.percentage);
	}
	else 
	{
		res.MapToInterval(0, 1);
		len = res.GetPixelsCount();
		for (size_t index = 0; index < len; index++)
		{					
			float q = res.GetData()[index];
			if (q > info.threshold)
			{
				bestCorners->push_back(index);
			}		
		}
	}

		
	return res;
}

/// <summary>
/// Corner Harris function
/// </summary>
/// <param name="a"></param>
/// <param name="b"></param>
/// <param name="c"></param>
/// <param name="k"></param>
/// <returns></returns>
float ImageUtils::CornerHarrisFunction(float a, float b, float c, float k)
{
	float detM = a * b - c * c;
	float traceM = a + b;

	return std::abs(detM - k * traceM * traceM);
}

/// <summary>
/// Corner Shi-Tomasi function
/// computed as minimal value of eigenvalues
/// </summary>
/// <param name="a"></param>
/// <param name="b"></param>
/// <param name="c"></param>
/// <param name="k"></param>
/// <returns></returns>
float ImageUtils::CornerShiTomasiFunction(float a, float b, float c, float k)
{
	float tmp = (-a - b);
	float D = tmp * tmp - 4 * (a * b - c * c);
	if (D < 0)
	{
		return 0;
	}

	D = std::sqrt(D);

	float l1 = (-tmp + D) / 2.0f;
	float l2 = (-tmp - D) / 2.0f;

	return std::min(l1, l2);
}

/// <summary>
/// Find positions of max values in input image
/// If max value is found, neighborhood (neighborhoodRadius) 
/// around it is ignored
/// Return array of indexes with max values in input
/// </summary>
/// <param name="input"></param>
/// <param name="neighborhoodRadius"></param>
/// <param name="percents"></param>
/// <returns></returns>
std::vector<size_t> ImageUtils::FindMaxValues(const Image2d<float> & input, const NeighborhoodKernel & k,
	float percents)
{

	struct Candidate
	{
		float q;
		size_t index;
		
		Candidate(float q, size_t index) :
			q(q), index(index) {}
	};
	
	std::vector<Candidate> candidates;
	
	for (size_t i = 0; i < input.GetPixelsCount(); i++)
	{
		float v = input.GetData()[i];		
		if (v != 0)
		{
			candidates.emplace_back(v, i);
		}
	}


	std::sort(candidates.begin(), candidates.end(),
		[](const Candidate & a, const Candidate & b) -> bool {
			return a.q > b.q;
	});

	//around each point, we have neighborhoodSize x neighborhoodSize "dead" area
	//so max total number of points is candidates count div neighborhood
	size_t count = static_cast<size_t>((candidates.size() / (k.radius * k.radius)) * percents);

	//filter close values
	std::vector<bool> visited(input.GetData().size(), false);

	int yy_start, yy_stop;
	int xx_start, xx_stop;
	
	std::vector<size_t> bestValues;
	
	for (auto & c : candidates)
	{		
		if (visited[c.index])
		{
			continue;
		}

		bestValues.push_back(c.index);

		if (bestValues.size() >= count)
		{
			break;
		}

		visited[c.index] = true;

		int x, y;
		input.GetPositionFromIndex(c.index, x, y);

		ImageUtils::CalcKernelStartStop(y, input.GetHeight() - 1, k, yy_start, yy_stop);
		ImageUtils::CalcKernelStartStop(x, input.GetWidth() - 1, k, xx_start, xx_stop);

		for (int yy = yy_start; yy <= yy_stop; yy++)
		{
			for (int xx = xx_start; xx <= xx_stop; xx++)
			{
				size_t ni = size_t(xx) + size_t(yy) * size_t(input.GetWidth());
				visited[ni] = true;
			}
		}
	}

	return bestValues;
}

//==============================================================================================
// Image tilling
//==============================================================================================

/// <summary>
/// Create tiles from input image
/// Tile size is tileW x tileH
/// if allSameSize is true, last tile is filled with 0 if it exceed image size
/// if allSameSize is false, last tile size is different and tile is fully filled
/// </summary>
/// <param name="input"></param>
/// <param name="tileW"></param>
/// <param name="tileH"></param>
/// <param name="allSameSize"></param>
/// <returns></returns>
template <typename T>
std::vector<ImageUtils::TileInfo<T>> ImageUtils::CreateTiles(const Image2d<T> & input,
	int tileW, int tileH,
	bool allSameSize)
{		

	int tx = 0;
	int ty = 0;

	std::vector<ImageUtils::TileInfo<T>> res;

	int y = 0;
	for (; y < input.GetHeight() - tileH; y += tileH)
	{		
		tx = 0;
		int x = 0;
		for (; x < input.GetWidth() - tileW; x += tileW)
		{
			ImageUtils::TileInfo<T> ti;
			ti.x = tx;
			ti.y = ty;
			ti.img = Image2d<T>(tileW, tileH, input.GetPixelFormat());
						
			size_t index = 0;
			for (int yy = y; yy < y + tileH; yy++)
			{
				for (int xx = x; xx < x + tileW; xx++)
				{
					auto val = input.GetPixelStart(xx, yy);

					for (size_t c = 0; c < input.GetChannelsCount(); c++)
					{
						ti.img.SetValue(val[c], c, index);
					}
					index++;
				}				
			}
			
			res.emplace_back(ti);
			tx++;
		}


		//last tile in row can have different
		//width	
		
		int lastTileW = (allSameSize) ? tileW : (input.GetWidth() - x);
		int endOffset = (allSameSize) ? (tileW - (input.GetWidth() - x)) : 0;
		
		ImageUtils::TileInfo<T> ti;
		ti.x = tx;
		ti.y = ty;
		ti.img = Image2d<T>(lastTileW, tileH, input.GetPixelFormat());
		
		size_t index = 0;
		for (int yy = y; yy < y + tileH; yy++)
		{
			for (int xx = x; xx < input.GetWidth(); xx++)
			{
				auto val = input.GetPixelStart(xx, yy);

				for (size_t c = 0; c < input.GetChannelsCount(); c++)
				{
					ti.img.SetValue(val[c], c, index);
				}
				index++;
			}
			index += endOffset;
		}

		res.emplace_back(ti);	
		
		ty++;
	}
	

	//========================================================

	//last tile in column can have different
	//width

	tx = 0;
	int x = 0;
	for (; x < input.GetWidth() - tileW; x += tileW)
	{
		int lastTileH = (allSameSize) ? tileH : (input.GetHeight() - y);		
		ImageUtils::TileInfo<T> ti;
		ti.x = tx;
		ti.y = ty;
		ti.img = Image2d<T>(tileW, lastTileH, input.GetPixelFormat());

		size_t index = 0;
		for (int yy = y; yy < input.GetHeight(); yy++)
		{
			for (int xx = x; xx < x + tileW; xx++)
			{
				auto val = input.GetPixelStart(xx, yy);

				for (size_t c = 0; c < input.GetChannelsCount(); c++)
				{
					ti.img.SetValue(val[c], c, index);
				}
				index++;
			}
		}

		res.emplace_back(ti);
		tx++;
	}

	//last tile in column/row can have different
	//width / height

	int lastTileW = (allSameSize) ? tileW : (input.GetWidth() - x);
	int lastTileH = (allSameSize) ? tileH : (input.GetHeight() - y);
	int endOffset = (allSameSize) ? (tileW - (input.GetWidth() - x)) : 0;
		
	ImageUtils::TileInfo<T> ti;
	ti.x = tx;
	ti.y = ty;
	ti.img = Image2d<T>(lastTileW, lastTileH, input.GetPixelFormat());

	size_t index = 0;
	for (int yy = y; yy < input.GetHeight(); yy++)
	{
		for (int xx = x; xx < input.GetWidth(); xx++)
		{
			auto val = input.GetPixelStart(xx, yy);

			for (size_t c = 0; c < input.GetChannelsCount(); c++)
			{
				ti.img.SetValue(val[c], c, index);
			}
			index++;
		}
		index += endOffset;
	}

	res.emplace_back(ti);

	return res;
}

//=================================================================================================

template void ImageUtils::CalcDerivatives(const Image2d<float> & input, int x, int y,
	float * dx, float * dy, float * dxx, float * dyy, float * dxy);
template void ImageUtils::CalcDerivatives(const Image2d<uint8_t> & input, int x, int y,
	float * dx, float * dy, float * dxx, float * dyy, float * dxy);


template Image2d<float> ImageUtils::DetectCornersHarris(const Image2d<uint8_t> & input, const CornerDetection & info,
	std::vector<size_t> * bestCorners);
template Image2d<float> ImageUtils::DetectCornersHarris(const Image2d<float> & input, const CornerDetection & info,
	std::vector<size_t> * bestCorners);


template void ImageUtils::DrawLine(Image2d<float> & input, const float * value, int x0, int y0, int x1, int y1);
template void ImageUtils::DrawLine(Image2d<uint8_t> & input, const uint8_t * value, int x0, int y0, int x1, int y1);

template size_t ImageUtils::FloodFill(Image2d<uint8_t> & input, int x, int y, uint8_t newVal, uint8_t valDifThreshold, NeighborhoodKernel::Connectivity con);
template size_t ImageUtils::FloodFill(Image2d<float> & input, int x, int y, float newVal, float valDifThreshold, NeighborhoodKernel::Connectivity con);


template size_t ImageUtils::FloodFill(Image2d<uint8_t> & input, int x, int y, uint8_t newVal, std::function<bool(uint8_t, uint8_t, size_t)> sameAreaCallback, NeighborhoodKernel::Connectivity con);
template size_t ImageUtils::FloodFill(Image2d<float> & input, int x, int y, float newVal, std::function<bool(float, float, size_t)> sameAreaCallback, NeighborhoodKernel::Connectivity con);


template Image2d<float> ImageUtils::SDF(const Image2d<uint8_t> & input, uint8_t bgValue);
template Image2d<float> ImageUtils::SDF(const Image2d<float> & input, float bgValue);


template std::list<std::vector<ImageUtils::Pixel>> ImageUtils::Vectorize(const Image2d<float> & input, float bgValue, double minLength2, bool joinLines);
template std::list<std::vector<ImageUtils::Pixel>> ImageUtils::Vectorize(const Image2d<uint8_t> & input, uint8_t bgValue, double minLength2, bool joinLines);

template std::vector<ImageUtils::TileInfo<float>> ImageUtils::CreateTiles(const Image2d<float> & input, int tileW, int tileH, bool allSameSize);
template std::vector<ImageUtils::TileInfo<uint8_t>> ImageUtils::CreateTiles(const Image2d<uint8_t> & input, int tileW, int tileH, bool allSameSize);
