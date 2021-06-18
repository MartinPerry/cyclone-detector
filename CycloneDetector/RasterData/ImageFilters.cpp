#include "./ImageFilters.h"

#include <memory>
#include <thread>

#include "./Image2d.h"
#include "./ImageGausFilter.h"
#include "./ImageUtils.h"

#include "../Utils/Logger.h"

#ifdef ENABLE_SIMD
#	include <immintrin.h>
#	include "../SIMD/avx_utils.h"
#endif

//https://fiveko.com/tutorials/image-processing/median-filter/

//=================================================================================================
// Convolution filter parent class
//=================================================================================================

/// <summary>
/// Callback for update of result from filter
/// Default is no update, no value is changed
/// </summary>
/// <param name="res"></param>
/// <param name="c"></param>
void ConvolutionFilter::ResultUpdateNone(float & res)
{
}

/// <summary>
/// Callback for update of result from filter
/// This will square every value
/// </summary>
/// <param name="res"></param>
/// <param name="c"></param>
void ConvolutionFilter::ResultUpdateSqr(float & res)
{
	res *= res;
}

/// <summary>
/// Default parent ctor
/// </summary>
/// <param name="radius"></param>
/// <param name="borderMode"></param>
/// <param name="filterResultCallback"></param>
ConvolutionFilter::ConvolutionFilter(NeighborhoodKernel k, ImageUtils::BorderMode borderMode,
	void(*filterResultCallback)(float & res)) :
	n(k),
	borderMode(borderMode),
	filterResultCallback(filterResultCallback),
	activeChannel(0)
{
}

ConvolutionFilter2d * ConvolutionFilter::GetAsConvolutionFilter2d()
{
	return nullptr;
}

ConvolutionFilterSeparable * ConvolutionFilter::GetAsConvolutionFilterSeparable()
{
	return nullptr;
}

const  ConvolutionFilter2d * ConvolutionFilter::GetAsConvolutionFilter2d() const
{
	return nullptr;
}

const ConvolutionFilterSeparable * ConvolutionFilter::GetAsConvolutionFilterSeparable() const
{
	return nullptr;
}

/// <summary>
/// Get size of filter
/// </summary>
/// <returns></returns>
int ConvolutionFilter::GetSize() const
{
	return n.GetSize();
}

void ConvolutionFilter::SetActiveChannel(size_t channel) const
{
	this->activeChannel = channel;
}

//=================================================================================================
// Separable convolution filter
// https://theailearner.com/2019/05/24/first-order-derivative-kernels-for-edge-detection/
// https://en.wikipedia.org/wiki/Separable_filter
// Y^T * X = 2D version of filter (row vector)
//=================================================================================================

ConvolutionFilterSeparable ConvolutionFilterSeparable::CreateSobelX()
{
	return ConvolutionFilterSeparable({ -1, 0, 1 }, { 1, 2, 1 }, Order::X_Y);
}

ConvolutionFilterSeparable ConvolutionFilterSeparable::CreateSobelY()
{
	return ConvolutionFilterSeparable({ 1, 2, 1 }, { -1, 0, 1 }, Order::Y_X);
}

ConvolutionFilterSeparable ConvolutionFilterSeparable::CreateScharrX()
{
	return ConvolutionFilterSeparable({ -1, 0, 1 }, { 3, 10, 3 }, Order::X_Y);
}

ConvolutionFilterSeparable ConvolutionFilterSeparable::CreateScharrY()
{
	return ConvolutionFilterSeparable({ 3, 10, 3 }, { -1, 0, 1 }, Order::Y_X);
}

ConvolutionFilterSeparable::ConvolutionFilterSeparable() :
	ConvolutionFilter(
		NeighborhoodKernel::CreateFromSize(0),
		ImageUtils::BorderMode::CLAMP,
		&ConvolutionFilter::ResultUpdateNone),
	mulX(0),
	mulY(0),
	order(Order::Y_X)
{
}

ConvolutionFilterSeparable::ConvolutionFilterSeparable(const std::vector<float> & data) :
	ConvolutionFilter(
		NeighborhoodKernel::CreateFromSize(0),
		ImageUtils::BorderMode::CLAMP,
		&ConvolutionFilter::ResultUpdateNone),
	mulX(1),
	mulY(1),
	kernelDataX(data),
	kernelDataY(data),
	order(Order::Y_X) //order is not important, x and y kernels are same
{
	this->n.radius = static_cast<int>(data.size() / 2.0);
}

ConvolutionFilterSeparable::ConvolutionFilterSeparable(const std::vector<float> & dataX,
	const std::vector<float> & dataY,
	Order order) :
	ConvolutionFilter(
		NeighborhoodKernel::CreateFromSize(0),
		ImageUtils::BorderMode::CLAMP,
		&ConvolutionFilter::ResultUpdateNone),
	mulX(1),
	mulY(1),
	kernelDataX(dataX),
	kernelDataY(dataY),
	order(order)
{
	this->n.radius = static_cast<int>(dataX.size() / 2.0);
}

ConvolutionFilterSeparable * ConvolutionFilterSeparable::GetAsConvolutionFilterSeparable()
{
	return this;
}

const ConvolutionFilterSeparable * ConvolutionFilterSeparable::GetAsConvolutionFilterSeparable() const
{
	return this;
}

template <typename T>
float ConvolutionFilterSeparable::RunX(const Image2d<T> & input, int x, int y) const
{
	float tmp = 0;

	if (input.GetChannelsCount() == 1)
	{
		int kn = n.GetSize();
		const T * v = input.GetPixelStart(x - n.radius, y);
		for (int i = 0; i < kn; i++)
		{
			tmp += kernelDataX[i] * v[i];
		}
	}
	else
	{
		tmp = 0;
		int ki = 0;
		for (int xx = x - n.radius; xx <= x + n.radius; xx++)
		{
			float kernelVal = kernelDataX[ki];
			ki++;

			const T * v = input.GetPixelStart(xx, y);

			tmp += kernelVal * v[this->activeChannel];
		}
	}

	tmp *= mulX;
	return tmp;
}

template <typename T>
float ConvolutionFilterSeparable::RunY(const Image2d<T> & input, int x, int y) const
{
	float tmp = 0;
	int ki = 0;
	for (int yy = y - n.radius; yy <= y + n.radius; yy++)
	{
		float kernelVal = kernelDataY[ki];
		ki++;

		const T * v = input.GetPixelStart(x, yy);
		tmp += kernelVal * v[this->activeChannel];
	}

	tmp *= mulY;
	return tmp;
}


template <typename T>
float ConvolutionFilterSeparable::RunWithBorderX(const Image2d<T> & input, int x, int y) const
{
	int ki = 0;
	float tmp = 0;

	for (int xx = x - n.radius; xx <= x + n.radius; xx++)
	{
		float kernelVal = kernelDataX[ki];
		ki++;

		const T * v = input.GetPixelStartWithBorder(xx, y, borderMode);
		tmp += kernelVal * v[this->activeChannel];
	}

	tmp *= mulX;
	return tmp;
}

template <typename T>
float ConvolutionFilterSeparable::RunWithBorderY(const Image2d<T> & input, int x, int y) const
{
	int ki = 0;
	float tmp = 0;

	for (int yy = y - n.radius; yy <= y + n.radius; yy++)
	{
		float kernelVal = kernelDataY[ki];
		ki++;

		const T * v = input.GetPixelStartWithBorder(x, yy, borderMode);
		tmp += kernelVal * v[this->activeChannel];
	}

	tmp *= mulY;

	return tmp;
}



//=================================================================================================
// Classic 2D convolution filter
//=================================================================================================

ConvolutionFilter2d ConvolutionFilter2d::CreatePrevitX()
{
	return ConvolutionFilter2d({
		1, 0, -1,
		1, 0, -1,
		1, 0, -1
		});
}

ConvolutionFilter2d ConvolutionFilter2d::CreatePrevitY()
{
	return ConvolutionFilter2d({
		1, 1, 1,
		0, 0, 0,
		-1, -1, -1
		});
}

ConvolutionFilter2d ConvolutionFilter2d::CreateSobelX()
{
	return ConvolutionFilter2d({
		-1, 0, 1,
		-2, 0, 2,
		-1, 0, 1
		});
}

ConvolutionFilter2d ConvolutionFilter2d::CreateSobelY()
{
	return ConvolutionFilter2d({
		-1, -2, -1,
		0, 0, 0,
		1, 2, 1
		});
}

ConvolutionFilter2d ConvolutionFilter2d::CreateScharrX()
{
	return ConvolutionFilter2d({
		3, 0, -3,
		10, 0, -10,
		3, 0, -3
		});
}

ConvolutionFilter2d ConvolutionFilter2d::CreateScharrY()
{
	return ConvolutionFilter2d({
		3, 10, 3,
		0, 0, 0,
		-3, -10, -3
		});
}

ConvolutionFilter2d ConvolutionFilter2d::CreateSharpen()
{
	return ConvolutionFilter2d({
		0, -1, 0,
		-1, 5, -1,
		0, -1, 0
		});
}

/// <summary>
/// https://dsp.stackexchange.com/questions/10605/kernels-to-compute-second-order-derivative-of-digital-image
/// </summary>
/// <returns></returns>
ConvolutionFilter2d ConvolutionFilter2d::CreateDxx()
{
	return ConvolutionFilter2d({
		   1, -2, 1,
		   2, -4, 2,
		   1, -2, 1
		});
}

/// <summary>
/// https://dsp.stackexchange.com/questions/10605/kernels-to-compute-second-order-derivative-of-digital-image
/// </summary>
/// <returns></returns>
ConvolutionFilter2d ConvolutionFilter2d::CreateDyy()
{
	return ConvolutionFilter2d({
		   1,  2,  1,
		  -2, -4, -2,
		   1,  2,  1
		});
}

/// <summary>
/// https://dsp.stackexchange.com/questions/10605/kernels-to-compute-second-order-derivative-of-digital-image
/// </summary>
/// <returns></returns>
ConvolutionFilter2d ConvolutionFilter2d::CreateDxy()
{
	return ConvolutionFilter2d({
			1, 0, -1,
			0, 0,  0,
		   -1, 0,  1
		});
}

/// <summary>
/// Create kernel of square size - entire square is filled with val
/// </summary>
/// <param name="radius"></param>
/// <param name="val"></param>
/// <returns></returns>
ConvolutionFilter2d ConvolutionFilter2d::CreateSquare(int radius, float val)
{
	int s = radius * 2 + 1;
	std::vector<float> df = std::vector<float>(s * s, val);
	return ConvolutionFilter2d(df);
}

/// <summary>
/// Creake circular kernel - circle is filled with val
/// </summary>
/// <param name="radius"></param>
/// <param name="val"></param>
/// <returns></returns>
ConvolutionFilter2d ConvolutionFilter2d::CreateCircle(int radius, float val)
{
	int s = radius * 2 + 1;
	std::vector<float> df = std::vector<float>(s * s, 0);

	int index = 0;
	for (int y = -radius; y <= radius; y++)
	{
		for (int x = -radius; x <= radius; x++)
		{
			if (x * x + y * y <= radius * radius)
			{				
				df[index] = val;				
			}			
			index++;
		}		
	}
	

	return ConvolutionFilter2d(df);
}

ConvolutionFilter2d ConvolutionFilter2d::CreateReflected(const ConvolutionFilter2d & f)
{
	if (f.GetSize() % 2 == 0)
	{
		//size of filter is even - cannot reflect
		return f;
	}

	std::vector<float> kernelData(f.kernelData.size(), 0);

	for (int y = 0; y < f.n.radius; y++)
	{
		for (int x = 0; x < f.n.radius; x++)
		{
			size_t indexFrom = (x + f.n.radius) + (y + f.n.radius) * f.GetSize();
			size_t indexTo = (-x + f.n.radius) + (-y + f.n.radius) * f.GetSize();

			float val1 = f.kernelData[indexFrom];
			float val2 = f.kernelData[indexTo];

			kernelData[indexTo] = val1;
			kernelData[indexFrom] = val2;
		}
	}

	return ConvolutionFilter2d(kernelData);
}

/// <summary>
/// default ctor
/// </summary>
ConvolutionFilter2d::ConvolutionFilter2d() :
	ConvolutionFilter(
		NeighborhoodKernel::CreateFromSize(0),
		ImageUtils::BorderMode::CLAMP,
		&ConvolutionFilter::ResultUpdateNone),
	mul(1)
{
}

/// <summary>
/// ctor with prepared kernel data
/// Kernel radius is calculated from data
/// Data size (linearized matrix) N x N, where N must be odd
/// </summary>
/// <param name="data"></param>
ConvolutionFilter2d::ConvolutionFilter2d(const std::vector<float> & data) :
	ConvolutionFilter(
		NeighborhoodKernel::CreateFromSize(0),
		ImageUtils::BorderMode::CLAMP,
		&ConvolutionFilter::ResultUpdateNone),
	mul(1),
	kernelData(data)
{
	double s = std::sqrt(kernelData.size());
	this->n.radius = static_cast<int>(s / 2.0);
}

/// <summary>
/// Get pointer to this
/// </summary>
/// <returns></returns>
ConvolutionFilter2d * ConvolutionFilter2d::GetAsConvolutionFilter2d()
{
	return this;
}

/// <summary>
/// Get pointer to this
/// </summary>
/// <returns></returns>
const ConvolutionFilter2d * ConvolutionFilter2d::GetAsConvolutionFilter2d() const
{
	return this;
}

/// <summary>
/// Get pointer to kernel value at [x, y]
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
float * ConvolutionFilter2d::GetValueAt(int x, int y)
{
	return &this->kernelData[x + y * this->n.GetSize()];
}

/// <summary>
/// Get pointer to kernel value at [x, y]
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
const float * ConvolutionFilter2d::GetValueAt(int x, int y) const
{
	return &this->kernelData[x + y * this->n.GetSize()];
}

/// <summary>
/// Create negated version (multiply all elementy by -1)
/// </summary>
void ConvolutionFilter2d::Negate()
{
	this->Multiply(-1.0f);
}

/// <summary>
/// Multiply all elements by constant val
/// </summary>
/// <param name="val"></param>
void ConvolutionFilter2d::Multiply(float val)
{
	for (size_t i = 0; i < this->kernelData.size(); i++)
	{
		this->kernelData[i] *= val;
	}
}

/// <summary>
/// Run filter on input for pixel at [x, y]
/// It will iterate pixel neighborhood, cumulate values based on kernel and
/// return new value for center pixel [x, y]
/// </summary>
/// <param name="input"></param>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
template <typename T>
float ConvolutionFilter2d::Run(const Image2d<T> & input, int x, int y) const
{
	int ki = 0;
	float tmp = 0;

	for (int yy = y - n.radius; yy <= y + n.radius; yy++)
	{
		for (int xx = x - n.radius; xx <= x + n.radius; xx++)
		{
			float kernelVal = kernelData[ki];
			ki++;

			const T * v = input.GetPixelStart(xx, yy);
			tmp += kernelVal * v[this->activeChannel];
		}
	}

	tmp *= mul;

	filterResultCallback(tmp);

	return tmp;
}

/// <summary>
/// Run filter on input for pixel at [x, y] -> respect border settings
/// pixels that will over/underflo image will be treated according to borderMode
/// (this is slower, because it requires the condition testing for each pixel)
/// It will iterate pixel neighborhood, cumulate values based on kernel and
/// return new value for center pixel [x, y]
/// </summary>
/// <param name="input"></param>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
template <typename T>
float ConvolutionFilter2d::RunWithBorder(const Image2d<T> & input, int x, int y) const
{
	int ki = 0;
	float tmp = 0;

	for (int yy = y - n.radius; yy <= y + n.radius; yy++)
	{
		for (int xx = x - n.radius; xx <= x + n.radius; xx++)
		{
			float kernelVal = kernelData[ki];
			ki++;

			const T * v = input.GetPixelStartWithBorder(xx, yy, borderMode);
			tmp += kernelVal * v[this->activeChannel];
		}
	}

	tmp *= mul;

	filterResultCallback(tmp);

	return tmp;
}

//=================================================================================================
// Image filtering methods
//=================================================================================================

/// <summary>
/// Binarize image (convert it to 2 values only)
/// using simple thresshold for each channel (stored in an array thressholds)
/// Binary values are in array, where [0] / [1] correspond to lower / higher value
/// 
/// Output image must have same dimension and channels count as input
/// Output and input can be the same
/// 
/// Example:
/// ImageFilters::ThressholdBinary(in, out, {125, 255, 12}, {{0, 255}, {25, 255}, {50, 100}}
/// -> 3 channel image
/// Threesshold in 1. channel -> 125: result = (value < 125 ) ? 0 : 255
/// Threesshold in 2. channel -> 255: result = (value < 255 ) ? 25 : 255
/// Threesshold in 3. channel -> 12 : result = (value < 12  ) ? 50 : 100
/// 
/// </summary>
/// <param name="input"></param>
/// <param name="output"></param>
/// <param name="thressholds">array of thresshold for each channel</param>
/// <param name="binaryValues">array of binary values that are used for each channel</param>
template <typename T>
void ImageFilters<T>::ThressholdBinary(const Image2d<T> & input, Image2d<T> & output,
	const T * thressholds, const std::array<T, 2> * binaryValues)
{
	if ((input.GetWidth() != output.GetWidth()) ||
		(input.GetHeight() != output.GetHeight()) ||
		(input.GetChannelsCount() != output.GetChannelsCount()))
	{
		MY_LOG_ERROR("Input and output image dimensions differ");
		return;
	}

	if (&input == &output)
	{
		MY_LOG_ERROR("Input and output image are the same");
		return;
	}


	size_t len = input.GetPixelsCount();

	for (size_t i = 0; i < len; i++)
	{
		const T * in = input.GetPixelStart(i);
		T * out = output.GetPixelStart(i);

		for (size_t c = 0; c < input.GetChannelsCount(); c++)
		{
			out[c] = (in[c] < thressholds[c]) ? binaryValues[c][0] : binaryValues[c][1];
		}
	}
}


/// <summary>
/// Adaptive threshold
/// 
/// Threshold Condition:
/// if (input - gaus(input) > delta) ? background : foreground
/// c = channel index
/// binaryValues[c][0] = background
/// binaryValues[c][1] = foreground
/// 
/// source: OpenCV 
///  https://docs.opencv.org/4.2.0/d7/d1b/group__imgproc__misc.html#ga72b913f352e4a1b1b397736707afcde3
/// 
/// </summary>
/// <param name="input"></param>
/// <param name="output"></param>
/// <param name="gauss"></param>
/// <param name="binaryValues"></param>
/// <param name="delta"></param>
template <typename T>
void ImageFilters<T>::ThressholdAdaptive(const Image2d<T> & input, Image2d<T> & output,
	ImageGausFilter<T> & gauss, const std::array<T, 2> * binaryValues, double delta)
{
	if ((input.GetWidth() != output.GetWidth()) ||
		(input.GetHeight() != output.GetHeight()) ||
		(input.GetChannelsCount() != output.GetChannelsCount()))
	{
		MY_LOG_ERROR("Input and output image dimensions differ");
		return;
	}

	if (&input == &output)
	{
		MY_LOG_ERROR("Input and output image are the same");
		return;
	}

	gauss.SetInput(input);


	bool hasOutput = (gauss.GetOutput() != nullptr);

	std::shared_ptr<Image2d<T>> gaussImg = nullptr;
	if (hasOutput == false)
	{
		//gaussImg = input.CreateEmpty();		
		gaussImg = std::make_shared<Image2d<T>>(input.GetWidth(),
			input.GetHeight(),
			std::vector<T>(input.GetData().size()),
			input.GetPixelFormat());

		gauss.SetOutput(gaussImg.get());
	}
	gauss.RunSeparable();

	size_t len = input.GetPixelsCount();

	for (size_t i = 0; i < len; i++)
	{
		const T * inI = input.GetPixelStart(i);
		const T * inG = gauss.GetOutput()->GetPixelStart(i);
		T * out = output.GetPixelStart(i);

		for (size_t c = 0; c < input.GetChannelsCount(); c++)
		{
			out[c] = ((double(inI[c]) - double(inG[c])) > -delta) ? binaryValues[c][0] : binaryValues[c][1];
		}
	}

	if (hasOutput == false)
	{
		gauss.SetOutput(nullptr);
	}
}

/// <summary>
/// Detect edges using Sobel filter
/// Use default edgeCallback
/// (input / output in place support)
/// </summary>
/// <param name="input"></param>
/// <param name="output"></param>
template <typename T>
void ImageFilters<T>::EdgeDetection(const Image2d<T> & input, Image2d<T> & output)
{
	std::function<T(float, float)> edgeCallback =
		[&](float a, float b) -> T {
		a = (a < 0) ? -a : a;
		b = (b < 0) ? -b : b;

		return ImageUtils::clamp_cast<T>(255 - (a + b) * 0.8);
	};

	ImageFilters<T>::EdgeDetection(input, output, edgeCallback);
}

/// <summary>
/// Detect edges using Sobel filter
/// edgeCallback is used to convert sobel dX and dY to final image value
/// (input / output in place support)
/// </summary>
/// <param name="input"></param>
/// <param name="output"></param>
/// <param name="edgeCallback"></param>
template <typename T>
void ImageFilters<T>::EdgeDetection(const Image2d<T> & input, Image2d<T> & output,
	std::function<T(float, float)> edgeCallback)
{
	if ((input.GetWidth() != output.GetWidth()) ||
		(input.GetHeight() != output.GetHeight()) ||
		(input.GetChannelsCount() != output.GetChannelsCount()))
	{
		MY_LOG_ERROR("Input and output image dimensions differ");
		return;
	}

	size_t len = output.GetPixelsCount();

	Image2d<float> dx, dy;

	if (len > 1048576)
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

	if (input.GetChannelsCount() == 1)
	{
		for (size_t i = 0; i < len; i++)
		{
			const float * a = dx.GetPixelStart(i);
			const float * b = dy.GetPixelStart(i);

			T * out = output.GetPixelStart(i);

			out[0] = edgeCallback(*a, *b);
		}
	}
	else
	{
		for (size_t i = 0; i < len; i++)
		{
			const float * a = dx.GetPixelStart(i);
			const float * b = dy.GetPixelStart(i);

			for (size_t c = 0; c < input.GetChannelsCount(); c++)
			{
				T v = edgeCallback(a[c], b[c]);
				output.SetValue(v, c, i);
			}
		}
	}
}

/// <summary>
/// Median filter optimized for binary images
/// High speed optimization is only for kernel radius <= 7 
/// (kernel size is 15x15)
/// Otherwise the speed is similar to openCV version
/// Background = 0
/// Foreground = 255
/// </summary>
/// <param name="input"></param>
/// <param name="output"></param>
/// <param name="k"></param>
/// <param name="fgValue"></param>
/// <param name="borderMode"></param>
template <typename T>
void ImageFilters<T>::MedianBinaryImage(const Image2d<T> & input, Image2d<T> & output,
	NeighborhoodKernel k,
	ImageUtils::BorderMode borderMode)
{
	if ((input.GetWidth() != output.GetWidth()) ||
		(input.GetHeight() != output.GetHeight()) ||
		(input.GetChannelsCount() != output.GetChannelsCount()))
	{
		MY_LOG_ERROR("Input and output image dimensions differ");
		return;
	}

	/*
	//https://stackoverflow.com/questions/46390516/greater-than-function-in-c
	uint32_t u_isgt(uint32_t a, uint32_t b)
	{
		uint32_t diff = ~(~b + a);
		return ((diff ^ ((b ^ a) & (diff ^ b))) >> 31) & 1;
	}
	*/

	const T * in = input.GetData().data();
	T * out = output.GetData().data();

	//+1 we want index = sum1 - 1 >= ... => sum1 >= .... + 1
	//[0],[1],[2] => sum1 = 3 = 3 values, but index is 2
	uint32_t center = (((2 * k.radius + 1) * (2 * k.radius + 1)) / 2) + 1;

	if (borderMode != ImageUtils::BorderMode::ZERO)
	{
		std::function<void(int x, int y)> mc = [&](int x, int y) {
			int sum = 0;
			int count = 0;
			for (int yy = y - k.radius; yy <= y + k.radius; yy++)
			{
				for (int xx = x - k.radius; xx <= x + k.radius; xx++)
				{
					const T val = input.GetPixelStartWithBorder(xx, yy, borderMode)[0];

					if (val == 255) sum++;										
					count++;
				}
			}

			if (sum > count / 2)
			{
				out[x + y * input.GetWidth()] = 255;
			}
		};

		for (int y = 0; y < input.GetHeight(); y++)
		{
			for (int x = 0; x < k.radius; x++)
			{
				mc(x, y);
			}
		}

		for (int y = 0; y < input.GetHeight(); y++)
		{
			for (int x = input.GetWidth() - k.radius; x < input.GetWidth(); x++)
			{
				mc(x, y);
			}
		}

		for (int y = 0; y < k.radius; y++)
		{
			for (int x = 0; x < input.GetWidth(); x++)
			{
				mc(x, y);
			}
		}

		for (int y = input.GetHeight() - k.radius; y < input.GetHeight(); y++)
		{
			for (int x = 0; x < input.GetWidth(); x++)
			{
				mc(x, y);
			}
		}

	}

#ifdef ENABLE_SIMD
	if constexpr (std::is_same<T, uint8_t>::value)
	{
		if (k.radius <= 7)
		{
			//for radius < 7 this is highly optimized
			//and process 16 values at once

			//can be SIMDed with 8bits
			//maximal numbers in sum is 15 * 15

			uint32_t sum = 0;
			for (int y = k.radius; y < input.GetHeight() - k.radius; y++)
			{
				int endX = input.GetWidth() - k.radius;
				int startX = k.radius;

				int endXSimd = endX;
				int len = endX - startX;
				len = len - len % 16;
				endXSimd = startX + len;

				//process 16 values at once
				for (int x = startX; x < endXSimd; x += 16)
				{
					__m128i sum2 = _mm_setzero_si128();
					for (int yy = y - k.radius; yy <= y + k.radius; yy++)
					{
						size_t offset = yy * input.GetWidth();
						for (int xx = x - k.radius; xx <= x + k.radius; xx++)
						{
							__m128i val = _mm_loadu_si128((const __m128i*)(in + xx + offset));
							auto valMask = _mm_and_si128(val, _mm_set1_epi8(1));
							sum2 = _mm_add_epi8(sum2, valMask);
						}
					}

					auto vv = _my_mm_cmpge_epu8(sum2, _mm_set1_epi8(center));
					//vv is 0 or 255 (-1 at signed mode)

					_mm_storeu_si128((__m128i*)(out + x + y * input.GetWidth()), vv);
				}

				//process remaining values
				for (int x = endXSimd; x < endX; x++)
				{
					sum = 0;
					for (int yy = y - k.radius; yy <= y + k.radius; yy++)
					{
						size_t offset = yy * input.GetWidth();
						for (int xx = x - k.radius; xx <= x + k.radius; xx++)
						{
							const T val = in[xx + offset];
							//gives 0 for val = 0
							//and 1 otherwise
							sum += (val & 1);
						}
					}

					if (sum >= center)
					{
						out[x + y * input.GetWidth()] = 255;
					}
				}
			}

			return;
		}
	}
#endif

	//universal version
	//however, this is slower or on par with openCV

	uint32_t sum = 0;
	for (int y = k.radius; y < input.GetHeight() - k.radius; y++)
	{
		int endX = input.GetWidth() - k.radius;
		int startX = k.radius;

		for (int x = startX; x < endX; x++)
		{
			sum = 0;			
			for (int yy = y - k.radius; yy <= y + k.radius; yy++)
			{
				for (int xx = x - k.radius; xx <= x + k.radius; xx++)
				{
					const T val = in[xx + yy * input.GetWidth()];
					if constexpr (std::is_same<T, uint8_t>::value)
					{
						//gives 0 for val = 0
						//and 1 otherwise
						sum += (val & 1);
					}
					else
					{
						if (val == 255) sum++;
					}

				}
			}

			//out[x + y * input.GetWidth()] = fgValue * u_isgt(sum, center);

			if (sum >= center)
			{
				out[x + y * input.GetWidth()] = 255;
			}
		}
	}

}

/// <summary>
/// Run full 2D convolution filter on input image
/// and store result in output
/// Type of input and output can differ
/// </summary>
/// <param name="input"></param>
/// <param name="output"></param>
/// <param name="f"></param>
template <typename T>
template <typename V>
void ImageFilters<T>::RunConvolutionFilter(const Image2d<T> & input, Image2d<V> & output,
	const ConvolutionFilter2d & f)
{
	if ((input.GetWidth() != output.GetWidth()) ||
		(input.GetHeight() != output.GetHeight()) ||
		(input.GetChannelsCount() != output.GetChannelsCount()))
	{
		MY_LOG_ERROR("Input and output image dimensions differ");
		return;
	}

	if constexpr (std::is_same<T, V>::value)
	{
		if (&input == &output)
		{
			MY_LOG_ERROR("Input and output image are the same");
			return;
		}
	}


	for (size_t c = 0; c < input.GetChannelsCount(); c++)
	{
		f.SetActiveChannel(c);

		//process borders
		if (f.borderMode != ImageUtils::BorderMode::ZERO)
		{
			ConvolutionFilter::LoopHelper(input, output,
				0, input.GetHeight(),
				0, f.n.radius,
				f, &ConvolutionFilter2d::RunWithBorder);

			ConvolutionFilter::LoopHelper(input, output,
				0, input.GetHeight(),
				input.GetWidth() - f.n.radius, input.GetWidth(),
				f, &ConvolutionFilter2d::RunWithBorder);

			ConvolutionFilter::LoopHelper(input, output,
				0, f.n.radius,
				0, input.GetWidth(),
				f, &ConvolutionFilter2d::RunWithBorder);

			ConvolutionFilter::LoopHelper(input, output,
				input.GetHeight() - f.n.radius, input.GetHeight(),
				0, input.GetWidth(),
				f, &ConvolutionFilter2d::RunWithBorder);
		}


		//process inner
		ConvolutionFilter::LoopHelper(input, output,
			f.n.radius, input.GetHeight() - f.n.radius,
			f.n.radius, input.GetWidth() - f.n.radius,
			f, &ConvolutionFilter2d::Run);
	}
}

bool ConvolutionFilter::enableSimd = true;

/// <summary>
/// Run separable convolution filter on input image
/// and store result in output
/// Type of input and output can differ
/// (Note: this requires temporary image of same size as input)
/// </summary>
/// <param name="input"></param>
/// <param name="output"></param>
/// <param name="f"></param>
template <typename T>
template <typename V>
void ImageFilters<T>::RunConvolutionFilter(const Image2d<T> & input, Image2d<V> & output,
	const ConvolutionFilterSeparable & f)
{
	if ((input.GetWidth() != output.GetWidth()) ||
		(input.GetHeight() != output.GetHeight()) ||
		(input.GetChannelsCount() != output.GetChannelsCount()))
	{
		MY_LOG_ERROR("Input and output image dimensions differ");
		return;
	}

	if constexpr (std::is_same<T, V>::value)
	{
		if (&input == &output)
		{
			MY_LOG_ERROR("Input and output image are the same");
			return;
		}
	}

	Image2d<V> temp = input.template CreateEmpty<V>();

	for (size_t c = 0; c < input.GetChannelsCount(); c++)
	{
		f.SetActiveChannel(c);

		if (f.order == ConvolutionFilterSeparable::Order::Y_X)
		{

			if (f.borderMode != ImageUtils::BorderMode::ZERO)
			{
				ConvolutionFilter::LoopHelper(input, temp,
					0, input.GetHeight(),
					0, f.n.radius,
					f, &ConvolutionFilterSeparable::RunWithBorderY);

				ConvolutionFilter::LoopHelper(input, temp,
					0, input.GetHeight(),
					input.GetWidth() - f.n.radius, input.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderY);

				ConvolutionFilter::LoopHelper(input, temp,
					0, f.n.radius,
					0, input.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderY);

				ConvolutionFilter::LoopHelper(input, temp,
					input.GetHeight() - f.n.radius, input.GetHeight(),
					0, input.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderY);
			}

			if ((input.GetChannelsCount() == 1) && (ConvolutionFilter::enableSimd))
			{
				f.LoopHelperSingleChannelY(input, temp,
					f.n.radius, input.GetHeight() - f.n.radius,
					f.n.radius, input.GetWidth() - f.n.radius);
			}
			else
			{
				ConvolutionFilter::LoopHelper(input, temp,
					f.n.radius, input.GetHeight() - f.n.radius,
					f.n.radius, input.GetWidth() - f.n.radius,
					f, &ConvolutionFilterSeparable::RunY);
			}

			if (f.borderMode != ImageUtils::BorderMode::ZERO)
			{
				ConvolutionFilter::LoopHelper(temp, output,
					0, temp.GetHeight(),
					0, f.n.radius,
					f, &ConvolutionFilterSeparable::RunWithBorderX);

				ConvolutionFilter::LoopHelper(temp, output,
					0, temp.GetHeight(),
					temp.GetWidth() - f.n.radius, temp.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderX);

				ConvolutionFilter::LoopHelper(temp, output,
					0, f.n.radius,
					0, temp.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderX);

				ConvolutionFilter::LoopHelper(temp, output,
					temp.GetHeight() - f.n.radius, temp.GetHeight(),
					0, temp.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderX);
			}

			if ((input.GetChannelsCount() == 1) && (ConvolutionFilter::enableSimd))
			{
				f.LoopHelperSingleChannelX(temp, output,
					f.n.radius, input.GetHeight() - f.n.radius,
					f.n.radius, input.GetWidth() - f.n.radius);
			}
			else
			{
				ConvolutionFilter::LoopHelper(temp, output,
					f.n.radius, input.GetHeight() - f.n.radius,
					f.n.radius, input.GetWidth() - f.n.radius,
					f, &ConvolutionFilterSeparable::RunX);
			}
		}
		else
		{

			if (f.borderMode != ImageUtils::BorderMode::ZERO)
			{
				ConvolutionFilter::LoopHelper(input, temp,
					0, input.GetHeight(),
					0, f.n.radius,
					f, &ConvolutionFilterSeparable::RunWithBorderX);

				ConvolutionFilter::LoopHelper(input, temp,
					0, input.GetHeight(),
					input.GetWidth() - f.n.radius, input.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderX);

				ConvolutionFilter::LoopHelper(input, temp,
					0, f.n.radius,
					0, input.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderX);

				ConvolutionFilter::LoopHelper(input, temp,
					input.GetHeight() - f.n.radius, input.GetHeight(),
					0, input.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderX);
			}

			if ((input.GetChannelsCount() == 1) && (ConvolutionFilter::enableSimd))
			{
				f.LoopHelperSingleChannelX(input, temp,
					f.n.radius, input.GetHeight() - f.n.radius,
					f.n.radius, input.GetWidth() - f.n.radius);
			}
			else
			{
				ConvolutionFilter::LoopHelper(input, temp,
					f.n.radius, input.GetHeight() - f.n.radius,
					f.n.radius, input.GetWidth() - f.n.radius,
					f, &ConvolutionFilterSeparable::RunX);
			}

			if (f.borderMode != ImageUtils::BorderMode::ZERO)
			{
				ConvolutionFilter::LoopHelper(temp, output,
					0, temp.GetHeight(),
					0, f.n.radius,
					f, &ConvolutionFilterSeparable::RunWithBorderY);

				ConvolutionFilter::LoopHelper(temp, output,
					0, temp.GetHeight(),
					temp.GetWidth() - f.n.radius, temp.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderY);

				ConvolutionFilter::LoopHelper(temp, output,
					0, f.n.radius,
					0, temp.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderY);

				ConvolutionFilter::LoopHelper(temp, output,
					temp.GetHeight() - f.n.radius, temp.GetHeight(),
					0, temp.GetWidth(),
					f, &ConvolutionFilterSeparable::RunWithBorderY);
			}

			if ((input.GetChannelsCount() == 1) && (ConvolutionFilter::enableSimd))
			{
				f.LoopHelperSingleChannelY(temp, output,
					f.n.radius, input.GetHeight() - f.n.radius,
					f.n.radius, input.GetWidth() - f.n.radius);
			}
			else
			{
				ConvolutionFilter::LoopHelper(temp, output,
					f.n.radius, input.GetHeight() - f.n.radius,
					f.n.radius, input.GetWidth() - f.n.radius,
					f, &ConvolutionFilterSeparable::RunY);
			}
		}
	}
}

//=================================================================================================
// Template specialisations
//=================================================================================================


template class ImageFilters<float>;
template class ImageFilters<uint8_t>;

template void ImageFilters<float>::RunConvolutionFilter(const Image2d<float> & input, Image2d<float> & output, const ConvolutionFilter2d & f);
template void ImageFilters<float>::RunConvolutionFilter(const Image2d<float> & input, Image2d<uint8_t> & output, const ConvolutionFilter2d & f);
template void ImageFilters<uint8_t>::RunConvolutionFilter(const Image2d<uint8_t> & input, Image2d<float> & output, const ConvolutionFilter2d & f);
template void ImageFilters<uint8_t>::RunConvolutionFilter(const Image2d<uint8_t> & input, Image2d<uint8_t> & output, const ConvolutionFilter2d & f);

template void ImageFilters<float>::RunConvolutionFilter(const Image2d<float> & input, Image2d<float> & output, const ConvolutionFilterSeparable & f);
template void ImageFilters<float>::RunConvolutionFilter(const Image2d<float> & input, Image2d<uint8_t> & output, const ConvolutionFilterSeparable & f);
template void ImageFilters<uint8_t>::RunConvolutionFilter(const Image2d<uint8_t> & input, Image2d<float> & output, const ConvolutionFilterSeparable & f);
template void ImageFilters<uint8_t>::RunConvolutionFilter(const Image2d<uint8_t> & input, Image2d<uint8_t> & output, const ConvolutionFilterSeparable & f);


//=================================================================================================
// Separable Filters

template float ConvolutionFilterSeparable::RunY(const Image2d<float> & input, int x, int y) const;
template float ConvolutionFilterSeparable::RunY(const Image2d<uint8_t> & input, int x, int y) const;
template float ConvolutionFilterSeparable::RunX(const Image2d<float> & input, int x, int y) const;
template float ConvolutionFilterSeparable::RunX(const Image2d<uint8_t> & input, int x, int y) const;

template float ConvolutionFilterSeparable::RunWithBorderY(const Image2d<float> & input, int x, int y) const;
template float ConvolutionFilterSeparable::RunWithBorderY(const Image2d<uint8_t> & input, int x, int y) const;
template float ConvolutionFilterSeparable::RunWithBorderX(const Image2d<float> & input, int x, int y) const;
template float ConvolutionFilterSeparable::RunWithBorderX(const Image2d<uint8_t> & input, int x, int y) const;

//=================================================================================================
// 2D filters

template float ConvolutionFilter2d::Run(const Image2d<float> & input, int x, int y) const;
template float ConvolutionFilter2d::Run(const Image2d<uint8_t> & input, int x, int y) const;

template float ConvolutionFilter2d::RunWithBorder(const Image2d<float> & input, int x, int y) const;
template float ConvolutionFilter2d::RunWithBorder(const Image2d<uint8_t> & input, int x, int y) const;
