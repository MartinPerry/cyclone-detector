#include "./ImageResize.h"

#include <cmath>

#include "./Image2d.h"
#include "./ImageUtils.h"

/// <summary>
/// Nearest neighbor resize
/// </summary>
/// <param name="input"></param>
/// <param name="newW"></param>
/// <param name="newH"></param>
/// <returns></returns>
template <typename T>
Image2d<T> ImageResize<T>::ResizeNearestNeighbor(const Image2d<T> & input, const ImageDimension & newDim)
{
	if ((input.GetWidth() == newDim.w) && (input.GetHeight() == newDim.h))
	{
		//no resize needed
		return input;
	}

	size_t newSize = size_t(newDim.w) * size_t(newDim.h) * input.GetChannelsCount();

	std::vector<T> resizedData;
	resizedData.resize(newSize);

	double xRatio = input.GetWidth() / static_cast<double>(newDim.w);
	double yRatio = input.GetHeight() / static_cast<double>(newDim.h);

	int px, py;
	for (int y = 0; y < newDim.h; y++)
	{
		py = static_cast<int>(std::floor(y * yRatio));
		size_t yW = size_t(y) * size_t(newDim.w);

		for (int x = 0; x < newDim.w; x++)
		{
			px = static_cast<int>(std::floor(x * xRatio));
			const T * tmp = input.GetPixelStart(px, py);

			size_t newIndexStart = (x + yW) * input.GetChannelsCount();

			for (size_t c = 0; c < input.GetChannelsCount(); c++)
			{
				resizedData[newIndexStart + c] = tmp[c];
			}
		}
	}

	return Image2d<T>(newDim.w,
		newDim.h,
		std::move(resizedData),
		input.GetPixelFormat());
}

/// <summary>
/// Resize image with average value from neighborhood with Kernel
/// It is recomended for downsampling !!!
/// </summary>
/// <param name="input"></param>
/// <param name="newDim"></param>
/// <param name="k"></param>
/// <returns></returns>
template <typename T>
Image2d<T> ImageResize<T>::ResizeAverage(const Image2d<T> & input, const ImageDimension & newDim, 
	const NeighborhoodKernel & k)
{
	if ((input.GetWidth() == newDim.w) && (input.GetHeight() == newDim.h))
	{
		//no resize needed
		return input;
	}

	size_t newSize = size_t(newDim.w) * size_t(newDim.h) * input.GetChannelsCount();

	std::vector<T> resizedData;
	resizedData.resize(newSize);


	double xRatio = input.GetWidth() / static_cast<double>(newDim.w);
	double yRatio = input.GetHeight() / static_cast<double>(newDim.h);

	std::vector<double> q;
	q.resize(input.GetChannelsCount());	

	int yy_start, yy_stop;
	int xx_start, xx_stop;

	for (int y = 0; y < newDim.h; y++)
	{
		size_t yW = size_t(y) * size_t(newDim.w);

		double y0 = (y + 0.5) * yRatio;
		//size_t yy_start = static_cast<size_t>(std::max(std::floor(y0) - kernelSize + 1, 0.0));
		//size_t yy_stop = std::min(yy_start + (2 * kernelSize - 1), input.GetHeight());

		ImageUtils::CalcKernelStartStop(y0, input.GetHeight() - 1, k, yy_start, yy_stop);

		for (int x = 0; x < newDim.w; x++)
		{

			double x0 = (x + 0.5) * xRatio;
			//size_t xx_start = static_cast<size_t>(std::max(std::floor(x0) - kernelSize + 1, 0.0));
			//size_t xx_stop = std::min(xx_start + (2 * kernelSize - 1), input.GetWidth());
			ImageUtils::CalcKernelStartStop(x0, input.GetWidth() - 1, k, xx_start, xx_stop);

			//==============================================================================
			//interpolation kernel

			double weight = 0;

			for (size_t c = 0; c < input.GetChannelsCount(); c++)
			{
				q[c] = 0;
			}


			for (int yy = yy_start; yy <= yy_stop; yy++)
			{				
				for (int xx = xx_start; xx <= xx_stop; xx++)
				{					
					const T * tmp = input.GetPixelStart(xx, yy);
					for (size_t c = 0; c < input.GetChannelsCount(); c++)
					{
						q[c] += tmp[c];
					}

					weight++;
				}
			}

			//==============================================================================

			size_t newIndexStart = (x + yW) * input.GetChannelsCount();

			ImageResize<T>::WriteValue(newIndexStart, input.GetChannelsCount(),
				q, weight,
				resizedData);
			
		}
	}


	return Image2d<T>(newDim.w,
		newDim.h,
		std::move(resizedData),
		input.GetPixelFormat());
}


/// <summary>
/// /// <summary>
/// Hermite resize
/// Taken from:
/// https://github.com/viliusle/Hermite-resize
/// https://github.com/viliusle/Hermite-resize/blob/master/src/hermite.js
/// Note: Our alpha is same as other channels, but in github code is different
/// </summary>
/// <param name="input"></param>
/// <param name="newW"></param>
/// <param name="newH"></param>
/// <returns></returns>
template <typename T>
Image2d<T> ImageResize<T>::ResizeHermite(const Image2d<T> & input, const ImageDimension & newDim)
{
	if ((input.GetWidth() == newDim.w) && (input.GetHeight() == newDim.h))
	{
		//no resize needed
		return input;
	}

	size_t newSize = newDim.w * newDim.h * input.GetChannelsCount();

	std::vector<T> resizedData;
	resizedData.resize(newSize);

	double xRatio = input.GetWidth() / static_cast<double>(newDim.w);
	double yRatio = input.GetHeight() / static_cast<double>(newDim.h);

	double xRatioHalf = std::ceil(xRatio / 2.0);
	double yRatioHalf = std::ceil(yRatio / 2.0);

	std::vector<double> g;
	g.resize(input.GetChannelsCount());

	int yy_start, yy_stop;
	int xx_start, xx_stop;

	for (int y = 0; y < newDim.h; y++)
	{
		size_t yW = size_t(y) * size_t(newDim.w);

		double y0 = (y) * yRatio;

		yy_start = static_cast<int>(std::floor(y0));
		yy_stop = static_cast<int>(std::ceil(y0 + yRatio));
		yy_stop = std::min(yy_stop, input.GetHeight());

		
		for (int x = 0; x < newDim.w; x++)
		{			
			double x0 = (x) * xRatio;

			xx_start = static_cast<int>(std::floor(x0));
			xx_stop = static_cast<int>(std::ceil(x0 + xRatio));
			xx_stop = std::min(xx_stop, input.GetWidth());

			//==============================================================================
			//interpolation kernel
			
			double weights = 0;

			for (size_t c = 0; c < input.GetChannelsCount(); c++)
			{
				g[c] = 0;
			}
										
			for (int yy = yy_start; yy < yy_stop; yy++)
			{
				double dy = (y0 - (yy)) / yRatioHalf;  //no need for abs() -> dy * dy				
				double wy2 = dy * dy;

				for (int xx = xx_start; xx < xx_stop; xx++)
				{
					double dx = (x0 - (xx)) / xRatioHalf; //no need for abs() -> dx * dx
					double w2 = wy2 + dx * dx;
					if (w2 >= 1)
					{
						//pixel too far
						continue;
					}
					
					double w = std::sqrt(w2);

					//hermite filter
					double weight = 2 * w2 * w - 3 * w2 + 1;
					
					const T * tmp = input.GetPixelStart(xx, yy);
					for (size_t c = 0; c < input.GetChannelsCount(); c++)
					{
						g[c] += weight * tmp[c];
					}
					
					weights += weight;
				}
			}

			//==============================================================================

			size_t newIndexStart = (x + yW) * input.GetChannelsCount();

			ImageResize<T>::WriteValue(newIndexStart, input.GetChannelsCount(),
				g, weights,
				resizedData);
		}
	}

	return Image2d<T>(newDim.w,
		newDim.h,
		std::move(resizedData),
		input.GetPixelFormat());
}


/// <summary>
/// Resize image with Lanczos
/// </summary>
/// <param name="input"></param>
/// <param name="newDim"></param>
/// <param name="k"></param>
/// <returns></returns>
template <typename T>
Image2d<T> ImageResize<T>::ResizeLanczos(const Image2d<T> & input, const ImageDimension & newDim, 
	const NeighborhoodKernel & k)
{
	if ((input.GetWidth() == newDim.w) && (input.GetHeight() == newDim.h))
	{
		//no resize needed
		return input;
	}

	size_t newSize = size_t(newDim.w) * size_t(newDim.h) * input.GetChannelsCount();

	std::vector<T> resizedData;
	resizedData.resize(newSize);


	double xRatio = input.GetWidth() / static_cast<double>(newDim.w);
	double yRatio = input.GetHeight() / static_cast<double>(newDim.h);

	std::vector<double> q;
	q.resize(input.GetChannelsCount());
	
	int yy_start, yy_stop;
	int xx_start, xx_stop;

	for (int y = 0; y < newDim.h; y++)
	{
		size_t yW = size_t(y) * size_t(newDim.w);
		
		double y0 = (y + 0.5) * yRatio;
		//size_t yy_start = static_cast<size_t>(std::max(std::floor(y0) - kernelSize + 1, 0.0));		
		//size_t yy_stop = std::min(yy_start + (2 * kernelSize - 1), input.GetHeight());				
		ImageUtils::CalcKernelStartStop(y0, input.GetHeight() - 1, k, yy_start, yy_stop);
		
		for (int x = 0; x < newDim.w; x++)
		{

			double x0 = (x + 0.5) * xRatio;
			//size_t xx_start = static_cast<size_t>(std::max(std::floor(x0) - kernelSize + 1, 0.0));
			//size_t xx_stop = std::min(xx_start + (2 * kernelSize - 1), input.GetWidth());
			ImageUtils::CalcKernelStartStop(x0, input.GetWidth() - 1, k, xx_start, xx_stop);

			//==============================================================================
			//interpolation kernel

			double weights = 0;

			for (size_t c = 0; c < input.GetChannelsCount(); c++)
			{
				q[c] = 0;
			}

			
			for (int yy = yy_start; yy <= yy_stop; yy++)
			{
							
				double dy = (y0 - yy);
				double wY = ImageResize::LanczosKernel1d(dy, static_cast<double>(k.radius));


				for (int xx = xx_start; xx <= xx_stop; xx++)
				{	
					double dx = (x0 - xx);					
					double wX = ImageResize::LanczosKernel1d(dx, static_cast<double>(k.radius));
					
					double wXY = wX * wY;

					const T * tmp = input.GetPixelStart(xx, yy);
					for (size_t c = 0; c < input.GetChannelsCount(); c++)
					{
						q[c] += tmp[c] * wXY;
					}

					weights += wXY;
				}												
			}

			//==============================================================================

			size_t newIndexStart = (x + yW) * input.GetChannelsCount();

			ImageResize<T>::WriteValue(newIndexStart, input.GetChannelsCount(),
				q, weights, 
				resizedData);
		}
	}
	

	return Image2d<T>(newDim.w,
		newDim.h,
		std::move(resizedData),
		input.GetPixelFormat());
}

//=========================================================================================

/// <summary>
/// Calculate 1D Lanczos filter with kernel size n
/// </summary>
/// <param name="x"></param>
/// <param name="n"></param>
/// <returns></returns>
template <typename T>
double ImageResize<T>::LanczosKernel1d(double x, double n)
{	
	static const double PI = 3.1415927410125732421875;
		
	if ((x > -0.00001) && (x < 0.00001)) return 1;
	
	if (std::abs(x) >= n) return 0;

	//normalized sinc
	//double sinc = std::sin(PI * x) / (PI * x);

	//combined equation
	//sinc * sin(PI * x/n) / (PI * x / n)

	double pix = PI * x;
	return n * (std::sin(pix) * std::sin(pix / n)) / (pix * pix);
}

//=========================================================================================
// Helper methods
//=========================================================================================


/// <summary>
/// Helper method to write result to resizedData
/// </summary>
/// <param name="newIndexStart"></param>
/// <param name="chanCount"></param>
/// <param name="vals"></param>
/// <param name="weight"></param>
/// <param name="resizedData"></param>
template <typename T>
void ImageResize<T>::WriteValue(size_t newIndexStart, size_t chanCount, 
	const std::vector<double> & vals, double weight,
	std::vector<T> & resizedData)
{
	for (size_t c = 0; c < chanCount; c++)
	{
		resizedData[newIndexStart + c] = ImageUtils::clamp_cast<T>((vals[c] / weight) + 0.5);		
	}
}

//=================================================================================================
// Template specialisations
//=================================================================================================

template class ImageResize<uint8_t>;
template class ImageResize<float>;

