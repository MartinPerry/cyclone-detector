#ifndef IMAGE_GAUS_FILTER_H
#define IMAGE_GAUS_FILTER_H

template <typename T>
class Image2d;

#include <array>
#include <vector>

#include "./ImageUtils.h"

template <typename T>
class ImageGausFilter
{
public:
	
	ImageGausFilter();
	ImageGausFilter(int radius, double weight = 1.0, double strength = 1.0);
	~ImageGausFilter();

	int GetRadius() const;
	void SetBorderMode(ImageUtils::BorderMode mode);

	void SetInput(const Image2d<T> & input);
	void SetOutput(Image2d<T> & output);
	void SetInput(const Image2d<T> * input);
	void SetOutput(Image2d<T> * output);
	Image2d<T> * GetOutput() const;

	void RunBruteForce();
	void RunSeparable();
	
	void RunBruteForce(int startX, int startY, size_t w, size_t h);
	void RunSeparable(int startX, int startY, size_t w, size_t h);

private:
	
	double strength;
	double weight;
	int radius;
	ImageUtils::BorderMode mode;
	
	const Image2d<T> * input;
	Image2d<T> * output;

	std::vector<float> kernel;

	void Init2DKernel();
	void Init1DKernel();
		
	
};


#endif