#ifndef IMAGE_RESIZE_H
#define IMAGE_RESIZE_H

template <typename T>
class Image2d;

#include <vector>

#include "./ImageUtils.h"

template <typename T>
class ImageResize 
{
public:
	
	static Image2d<T> ResizeNearestNeighbor(const Image2d<T> & input, const ImageDimension & newDim);
	
	static Image2d<T> ResizeAverage(const Image2d<T> & input, const ImageDimension & newDim, 
		const NeighborhoodKernel & k);
	
	static Image2d<T> ResizeHermite(const Image2d<T> & input, const ImageDimension & newDim);
	
	static Image2d<T> ResizeLanczos(const Image2d<T> & input, const ImageDimension & newDim, 
		const NeighborhoodKernel & k);

protected:
				
	static void WriteValue(size_t newIndexStart, size_t chanCount, 
		const std::vector<double> & vals, double weight,
		std::vector<T> & resizedData);

	static double LanczosKernel1d(double x, double n);

};

#endif
