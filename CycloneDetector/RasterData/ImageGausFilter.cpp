#include "./ImageGausFilter.h"

#include <cmath>

#include "./Image2d.h"

#include "./ImageFilters.h"

#include "../Macros.h"

//=================================================================================================

/// <summary>
/// default ctor
/// weight = 1
/// strength = 1
/// radius = 1
/// border mode = ENLARGE
/// </summary>
template <typename T>
ImageGausFilter<T>::ImageGausFilter() : 
	strength(1.0),
	weight(1.0),	
	radius(1),
	mode(ImageUtils::BorderMode::ENLARGE),
	input(nullptr),
	output(nullptr)
{	
}

/// <summary>
/// ctor
/// border mode = ENLARGE
/// </summary>
/// <param name="radius"></param>
/// <param name="weight"></param>
/// <param name="strength"></param>
template <typename T>
ImageGausFilter<T>::ImageGausFilter(int radius, double weight, double strength) : 
	strength(strength),
	weight(weight),	
	radius(radius),
	mode(ImageUtils::BorderMode::ENLARGE),	
	input(nullptr),
	output(nullptr)
{	
}

/// <summary>
/// dtor
/// </summary>
template <typename T>
ImageGausFilter<T>::~ImageGausFilter()
{
	
}

//=================================================================================================

/// <summary>
/// Get filter radius
/// </summary>
template <typename T>
int ImageGausFilter<T>::GetRadius() const
{
	return this->radius;
}

/// <summary>
/// Set mode for treating image borders
/// </summary>
/// <param name="mode"></param>
template <typename T>
void ImageGausFilter<T>::SetBorderMode(ImageUtils::BorderMode mode)
{
	this->mode = mode;
}


template <typename T>
void ImageGausFilter<T>::SetInput(const Image2d<T> & input)
{
	this->input = &input;
}

template <typename T>
void ImageGausFilter<T>::SetOutput(Image2d<T> & output)
{
	this->output = &output;
}

template <typename T>
void ImageGausFilter<T>::SetInput(const Image2d<T> * input)
{
	this->input = input;
}

template <typename T>
void ImageGausFilter<T>::SetOutput(Image2d<T> * output)
{
	this->output = output;
}

template <typename T>
Image2d<T> * ImageGausFilter<T>::GetOutput() const
{
	return this->output;
}

/// <summary>
/// Init kernel for brute-force version of filter
/// </summary>
template <typename T>
void ImageGausFilter<T>::Init2DKernel()
{
	
	int maxSize = (this->radius * 2 + 1) * (this->radius * 2 + 1);
	this->kernel.clear();
	this->kernel.resize(maxSize);
		
	//Full 2D filter
	//http://code.msdn.microsoft.com/Calculating-Gaussian-c7c95d2c
		
	//1D filter
	//http://en.wikipedia.org/wiki/Gaussian_blur

	   
	float sumTotal = 0;       
    double distance = 0;  
   
	double kernelMul = 1.0 / (2.0 * 3.14 * weight * weight);  
	
	int index = 0;
	for (int filterY = -this->radius;  filterY <= this->radius; filterY++)  
    { 
        for (int filterX = -this->radius;  filterX <= this->radius; filterX++)  
        { 
            distance = filterX * (1.0 - this->strength) * filterX * (1.0 - this->strength);
			distance += filterY * (1.0 - this->strength) * filterY * (1.0 - this->strength);
			distance /= (2.0 * (this->weight * this->weight));  
   
            this->kernel[index] =  static_cast<float>(kernelMul * std::exp(-distance));  
   
            sumTotal += this->kernel[index];  
			index++;
        }  
    }  
  	
    for (int i = 0; i < maxSize; i++)  
    {  
        this->kernel[i] /= sumTotal;  		
    }        
}

/// <summary>
/// Init kernel for separable version of filter
/// </summary>
template <typename T>
void ImageGausFilter<T>::Init1DKernel()
{
	
	int maxSize = this->radius * 2 + 1;
	this->kernel.clear();
	this->kernel.resize(maxSize);
		
		
	//1D filter
	//http://en.wikipedia.org/wiki/Gaussian_blur
	 
   
	float sumTotal = 0;
     
    double distance = 0;     
	double kernelMul = 1.0 / std::sqrt(2.0 * 3.14 * this->weight * this->weight);  
	
	int index = 0;
	
	for (int filterX = -this->radius;  filterX <= this->radius; filterX++)  
    { 
		distance = (filterX * (1.0 - this->strength) * filterX * (1.0 - this->strength));
		distance /= (2.0 * (this->weight * this->weight));  
   
        this->kernel[index] = static_cast<float>(kernelMul * exp(-distance));
   
        sumTotal += this->kernel[index];  
		index++;        
    }  
  	
    for (int i = 0; i < maxSize; i++)  
    {  
        this->kernel[i] /= sumTotal;  		
    }         
}


//=================================================================================================


template <typename T>
void ImageGausFilter<T>::RunBruteForce()
{
	this->RunBruteForce(0, 0, input->GetWidth(), input->GetHeight());	
}

template <typename T>
void ImageGausFilter<T>::RunBruteForce(int startX, int startY, 
	size_t w, size_t h)
{
	if ((input == nullptr) || (output == nullptr))
	{		
		return;
	}

	this->Init2DKernel();
	
	Image2d<T> * tmp = output;
	if (input == tmp)
	{
		tmp = new Image2d<T>(
					input->GetWidth(),
					input->GetHeight(),
					input->GetData(),
					input->GetPixelFormat()
				);		
	}
	

	ConvolutionFilter2d gaus2d(this->kernel);
	gaus2d.borderMode = this->mode;

	ImageFilters<T>::RunConvolutionFilter(*input, *tmp, gaus2d);
	

	if (input == output)
	{
		*output = std::move(*tmp);
		SAFE_DELETE(tmp);
	}
	
}


//=================================================================================================

/// <summary>
/// Run separable Gaus filter on entire image
/// if input == output, input data are overwritten,
/// but helper temporary image is created during computation
/// 
/// </summary>
/// <param name="input"></param>
/// <param name="output"></param>
template <typename T>
void ImageGausFilter<T>::RunSeparable()
{
	this->RunSeparable(0, 0, input->GetWidth(), input->GetHeight());	
	
}

/// <summary>
/// Run separable Gaus filter on sub-simage of input
/// Subimage is given by its top-left corner [startX, startY]
/// and size w x h
/// If we want to run entire image, use start as [0, 0] and image w, h
/// 
/// if input == output, input data are overwritten,
/// but helper temporary image is created during computation
/// </summary>
/// <param name="input"></param>
/// <param name="output"></param>
/// <param name="startX"></param>
/// <param name="startY"></param>
/// <param name="w"></param>
/// <param name="h"></param>
template <typename T>
void ImageGausFilter<T>::RunSeparable(int startX, int startY, 
	size_t w, size_t h)
{
	if ((input == nullptr) || (output == nullptr))
	{
		return;
	}

	this->Init1DKernel();
	
	Image2d<T> * tmp = output;
	if (input == tmp)
	{
		tmp = new Image2d<T>(
				input->GetWidth(),
				input->GetHeight(),
				input->GetData(),
				input->GetPixelFormat()
			);
	}


	//http://www.programming-techniques.com/2013/03/gaussian-blurring-using-separable.html
	//https://fiveko.com/tutorials/image-processing/gaussian-blur-filter/


	ConvolutionFilterSeparable gausSep(this->kernel);
	gausSep.borderMode = this->mode;

	ImageFilters<T>::RunConvolutionFilter(*input, *tmp, gausSep);
		
	if (input == output)
	{
		*output = std::move(*tmp);
		SAFE_DELETE(tmp);
	}	
}


//=================================================================================================

template class ImageGausFilter<uint8_t>;
template class ImageGausFilter<float>;
