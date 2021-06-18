#include "./Image2d.h"

#include <algorithm>

#include "../Compression/PNGLoader.h"
#include "../Compression/3rdParty/lodepng.h"

#include "../Utils/Logger.h"
#include "../Utils/FileWrapper.h"

#include "../Math/MathUtils.h"

#include "./ImageGausFilter.h"

#ifdef HAVE_OPENCV
#	include <opencv2/highgui.hpp>
#	ifdef _DEBUG
#		pragma comment(lib, "opencv_world452d.lib")
#	else
#		pragma comment(lib, "opencv_world452.lib")
#	endif
#endif

#ifdef ENABLE_SIMD
#	include <immintrin.h>
#	include "../SIMD/avx_utils.h"
#endif

#ifdef _MSC_VER
#	ifndef my_fopen 
#		define my_fopen(a, b, c) fopen_s(a, b, c)	
#	endif	
#else
#	ifndef my_fopen 
#		define my_fopen(a, b, c) (*a = fopen(b, c))
#	endif	
#endif

#ifdef MYMATH_NAMESPACE
using namespace MYMATH_NAMESPACE;
#endif

//=================================================================================================
// static factories
//=================================================================================================

/// <summary>
/// Load data from RAW file
/// User must specify image dimension and format,
/// loaded file is loaded directly as it is
/// </summary>
/// <param name="w"></param>
/// <param name="h"></param>
/// <param name="pf"></param>
/// <param name="fileName"></param>
/// <returns></returns>
template <typename T>
Image2d<T> Image2d<T>::CreateFromRawFile(int w, int h, ColorSpace::PixelFormat pf, const char * fileName)
{
	Image2d<T> img(w, h, pf);

	RawFile f = RawFile(fileName, "rb");
	if (f.GetSize() != img.data.size() * sizeof(T))
	{
		return img;
	}
	f.Read(&img.data[0], sizeof(T), img.data.size());
	return img;
}


//=================================================================================================
// ctors & dtor
//=================================================================================================

/// <summary>
/// Create empty image with size 0
/// </summary>
template <typename T>
Image2d<T>::Image2d() : 
	dim({ 0, 0 }),
	pf(ColorSpace::PixelFormat::NONE),
	channelsCount(ColorSpace::GetChannelsCount(ColorSpace::PixelFormat::NONE))
{
}

/// <summary>
/// Create image form file with fileName
/// Only supported compression of file are PNG and JPG
/// (Files are loaded as they are. If they contains color profile,
/// it is ignored)
/// </summary>
/// <param name="fileName"></param>
template <typename T>
Image2d<T>::Image2d(const char * fileName) :
	dim({ 0, 0 }),
	pf(ColorSpace::PixelFormat::NONE),
	channelsCount(ColorSpace::GetChannelsCount(ColorSpace::PixelFormat::NONE))
{
	RawFile f = RawFile(fileName, "rb");
	if (f.GetRawFilePtr() == nullptr)
	{
		MY_LOG_ERROR("File %s not found", fileName);
		return;
	}

	uint8_t header[8];
	f.Read(header, sizeof(uint8_t), 8);
	f.Seek(0, SEEK_SET);

	if ((header[0] == 137) && (header[1] == 'P')) // N G \r \n \032 \n
	{
		//PNG
		PNGLoader png;
		png.SetKeepPalette(false); //we do not want palette - just a simple load
		auto dPng = png.DecompressFromFile(&f);

		this->dim.w = dPng.w;
		this->dim.h = dPng.h;
		
		switch (dPng.channelsCount)
		{
		case 1:
			this->pf = ColorSpace::PixelFormat::GRAY;
			break;
		case 3:
			this->pf = ColorSpace::PixelFormat::RGB;
			break;
		case 4:
			this->pf = ColorSpace::PixelFormat::RGBA;
			break;
		default:
			this->pf = ColorSpace::PixelFormat::NONE;
			break;
		}

		this->channelsCount = dPng.channelsCount;
		if constexpr (std::is_same<T, uint8_t>::value)
		{
			this->data = std::move(dPng.data);
		}
		else 
		{
			this->data.resize(dPng.data.size());

			for (size_t i = 0; i < this->data.size(); i++)
			{
				this->data[i] = static_cast<T>(dPng.data[i] / 255.0);
			}
		}
	}
		
}


/// <summary>
/// Create empty image with size w x h and given pixel format
/// </summary>
/// <param name="w"></param>
/// <param name="h"></param>
/// <param name="pf"></param>
template <typename T>
Image2d<T>::Image2d(int w, int h, ColorSpace::PixelFormat pf) :
	dim({ w, h }),	
	data(std::vector<T>(w * h * ColorSpace::GetChannelsCount(pf), 0)),
	pf(pf),
	channelsCount(ColorSpace::GetChannelsCount(pf))
{		
}

/// <summary>
/// Create image with size w x h, filled with data and given pixel format
/// </summary>
/// <param name="w"></param>
/// <param name="h"></param>
/// <param name="data"></param>
/// <param name="pf"></param>
template <typename T>
Image2d<T>::Image2d(int w, int h, const std::vector<T> & data, ColorSpace::PixelFormat pf) :
	dim({ w, h }),
	data(data),
	pf(pf),
	channelsCount(ColorSpace::GetChannelsCount(pf))
{
}

/// <summary>
/// Create image with size w x h, filled with data and given pixel format
/// </summary>
/// <param name="w"></param>
/// <param name="h"></param>
/// <param name="data"></param>
/// <param name="pf"></param>
template <typename T>
Image2d<T>::Image2d(int w, int h, const T * rawData, ColorSpace::PixelFormat pf) :
	dim({ w, h }),
	data(rawData, rawData + w * h * ColorSpace::GetChannelsCount(pf)),
	pf(pf),
	channelsCount(ColorSpace::GetChannelsCount(pf))
{
}

template <typename T>
Image2d<T>::Image2d(int w, int h, const T ** rawData, ColorSpace::PixelFormat pf) :
	dim({ w, h }),	
	pf(pf),
	channelsCount(ColorSpace::GetChannelsCount(pf))
{
	for (int i = 0; i < h; i++)
	{
		data.insert(data.end(), rawData[i], rawData[i] + w * channelsCount);
	}
}

/// <summary>
/// Create image with size w x h, filled with data and given pixel format
/// data are moved
/// </summary>
/// <param name="w"></param>
/// <param name="h"></param>
/// <param name="data"></param>
/// <param name="pf"></param>
template <typename T>
Image2d<T>::Image2d(int w, int h, std::vector<T> && data, ColorSpace::PixelFormat pf) :
	dim({ w, h }),	
	data(std::move(data)),
	pf(pf),
	channelsCount(ColorSpace::GetChannelsCount(pf))
{	
}

template <typename T>
Image2d<T>::Image2d(const Image2d<T> & other) : 
	dim(other.dim),	
	data(other.data),
	pf(other.pf),
	channelsCount(other.channelsCount)
{
}

/// <summary>
/// move ctor
/// </summary>
/// <param name="other"></param>
template <typename T>
Image2d<T>::Image2d(Image2d<T> && other) noexcept :
	dim(other.dim),
	data(std::move(other.data)),
	pf(other.pf),
	channelsCount(other.channelsCount)
{
	other.Release();
}


#ifdef HAVE_OPENCV

/// <summary>
/// Copy data from cv::Mat to Image
/// Data are treated as GRAY/RGB/RGBA
/// Image owns the data
/// 
/// Source: https://stackoverflow.com/questions/26681713/convert-mat-to-array-vector-in-opencv
/// </summary>
/// <param name="cvMat"></param>
template <typename T>
Image2d<T>::Image2d(const cv::Mat & cvMat) :
	dim({ cvMat.cols, cvMat.rows }),	
	channelsCount(cvMat.channels())
{	
	if (channelsCount == 1) pf = ColorSpace::PixelFormat::GRAY;
	else if (channelsCount == 3) pf = ColorSpace::PixelFormat::RGB;
	else if (channelsCount == 4) pf = ColorSpace::PixelFormat::RGBA;
	
	if (cvMat.isContinuous())
	{		
		data.assign((T*)cvMat.data, (T*)cvMat.data + cvMat.total() * channelsCount);		
	}
	else 
	{
		for (int i = 0; i < cvMat.rows; i++)
		{
			data.insert(data.end(), cvMat.ptr<T>(i), cvMat.ptr<T>(i) + cvMat.cols);
		}
	}
}

#endif

/// <summary>
/// dtor
/// </summary>
template <typename T>
Image2d<T>::~Image2d()
{	
}

/// <summary>
/// Manually release data
/// dimensions of image are set to 0
/// </summary>
template <typename T>
void Image2d<T>::Release() noexcept
{
	this->data.clear();
	this->channelsCount = 0;
	this->dim.w = 0;
	this->dim.h = 0;
	this->pf = ColorSpace::PixelFormat::NONE;
}

//=================================================================================================
// Operators section
//=================================================================================================


template <typename T>
Image2d<T> & Image2d<T>::operator=(const Image2d<T> & other)
{
	this->dim = other.dim;	
	this->data = other.data;
	this->pf = other.pf;
	this->channelsCount = other.channelsCount;
	return *this;
}

template <typename T>
Image2d<T> & Image2d<T>::operator=(Image2d<T> && other) noexcept
{
	this->dim = other.dim;	
	this->data = std::move(other.data);	
	this->pf = other.pf;
	this->channelsCount = other.channelsCount;

	other.Release();

	return *this;
}

//=================================================================================================
// Creators section
// - create new image from existing
//=================================================================================================

/// <summary>
/// Create empty image with same parametrs as the current one
/// </summary>
/// <returns></returns>
template <typename T>
template <typename V>
Image2d<V> Image2d<T>::CreateEmpty() const
{	
	return Image2d<V>(this->GetWidth(),
		this->GetHeight(),
		std::vector<V>(this->data.size(), 0),
		this->pf);
}

/// <summary>
/// Create deep copy of current image
/// </summary>
/// <returns></returns>
template <typename T>
Image2d<T> Image2d<T>::CreateDeepCopy() const
{
	return Image2d<T>(this->GetWidth(),
		this->GetHeight(),
		this->data,
		this->pf);
}

/// <summary>
/// Create new image using only single channel of the current one
/// Output image is PixelFormat::GRAY
/// </summary>
/// <param name="channelIndex"></param>
/// <returns></returns>
template <typename T>
Image2d<T> Image2d<T>::CreateFromChannel(size_t channelIndex) const
{
	size_t len = size_t(this->dim.w) * size_t(this->dim.h);

	std::vector<T> channelData;
	channelData.resize(len);
	
	for (size_t i = 0; i < len; i++)
	{
		const T * tmp = this->GetPixelStart(i);

		channelData[i] = tmp[channelIndex];
	}

	return Image2d<T>(this->GetWidth(),
		this->GetHeight(),
		std::move(channelData),
		ColorSpace::PixelFormat::GRAY);
}

/// <summary>
/// Create new image from the current image sub-image
/// sub-image starts at [x, y] and has given size 
/// No checks are performed - x, y, and size must be in currents
/// image bounds
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="size"></param>
/// <returns></returns>
template <typename T>
Image2d<T> Image2d<T>::CreateSubImage(int x, int y, const ImageDimension & size) const
{
	size_t len = size_t(size.w) * size_t(size.h) * this->channelsCount;

	std::vector<T> subData;
	subData.resize(len);

	size_t index = 0;
	for (int yy = y; yy < y + size.h; yy++)
	{
		for (int xx = x; xx < x + size.w; xx++)
		{
			const T * tmp = this->GetPixelStart(xx, yy);
			for (size_t c = 0; c < this->channelsCount; c++)
			{
				subData[index] = tmp[c];
				index++;
			}
		}
	}

	return Image2d<T>(size.w,
		size.h,
		std::move(subData),
		this->pf);
}

/// <summary>
/// Create new Image2d from current
/// but cast data from T to V
/// 
/// Use to cast image from:
/// float -> uint8_t
/// uint8_t -> float
/// 
/// Data must be in correct range
/// </summary>
/// <returns></returns>
template <typename T>
template <typename V>
Image2d<V> Image2d<T>::CreateAs() const
{	
	std::vector<V> d;
	d.resize(this->data.size());

	size_t dataSize4 = this->data.size() - this->data.size() % 4;

	for (size_t i = 0; i < dataSize4; i += 4)
	{		
		d[i] = static_cast<V>(this->data[i]);
		d[i + 1] = static_cast<V>(this->data[i + 1]);
		d[i + 2] = static_cast<V>(this->data[i + 2]);
		d[i + 3] = static_cast<V>(this->data[i + 3]);
	}

	for (size_t i = dataSize4; i < this->data.size(); i++)
	{
		d[i] = static_cast<V>(this->data[i]);
	}

	return Image2d<V>(this->GetWidth(),
		this->GetHeight(),
		std::move(d),
		this->pf);
}

/// <summary>
/// Create new Image2d from current
/// but cast data from T to V
/// and map them to interval [start, end] for each channel
/// For mapping, min and max of each channel is computed
/// 
/// Use to cast image from:
/// float -> uint8_t
/// uint8_t -> float
/// 
/// </summary>
/// <param name="start"></param>
/// <param name="end"></param>
/// <returns></returns>
template <typename T>
template <typename V>
Image2d<V> Image2d<T>::CreateAsMapped(V start, V end) const
{
	std::vector<V> d;
	d.resize(this->data.size());

	for (size_t c = 0; c < this->channelsCount; c++)
	{
		T min, max;
		this->FindMinMax(c, min, max);

		for (size_t i = c; i < this->data.size(); i += this->channelsCount)
		{
			float mappedVal = MathUtils::MapRange<float>(min, max, start, end, this->data[i]);

			d[i] = static_cast<V>(mappedVal);
		}
	}

	return Image2d<V>(this->GetWidth(),
		this->GetHeight(),
		std::move(d),
		this->pf);
}


//=================================================================================================
// Setters
//=================================================================================================

template <typename T>
void Image2d<T>::SetPixelFormat(ColorSpace::PixelFormat pf) noexcept
{
	this->pf = pf;
	this->channelsCount = ColorSpace::GetChannelsCount(pf);
}

/// <summary>
/// Save file to PNG or RAW
/// </summary>
/// <param name="fileName"></param>
template <typename T>
void Image2d<T>::Save(const char * fileName) const
{
	
	size_t len = strlen(fileName);

	if ((len > 4) &&
		(fileName[len - 4] == '.') && (fileName[len - 3] == 'r') &&
		(fileName[len - 2] == 'a') && (fileName[len - 1] == 'w'))
	{
		FILE * f = nullptr;
		my_fopen(&f, fileName, "wb");
		fwrite(this->data.data(), sizeof(T), this->data.size(), f);
		fclose(f);

		return;
	}


	uint8_t* rawData = nullptr;

	if constexpr (std::is_same<T, uint8_t>::value)
	{
		//PngSaver can modify data by removing premultiplied alpha
		//but we dont use this here
		rawData = (uint8_t *)(this->data.data());
	}
	else 
	{
		//data must be mapped to 0 - 1 interval
		//for floats -> we convert them from 0 - 1 to  0 - 255

		uint8_t * d = new uint8_t[this->data.size()];
		for (size_t i = 0; i < this->data.size(); i++)
		{
			d[i] = ImageUtils::clamp_cast<uint8_t>(this->data[i] * 255.0f);
		}

		rawData = d;
	}
	
	if ((len > 4) &&
		(fileName[len - 4] == '.') && (fileName[len - 3] == 'p') &&
		(fileName[len - 2] == 'n') && (fileName[len - 1] == 'g'))
	{				
		if (this->pf == ColorSpace::PixelFormat::GRAY)
		{
			lodepng::encode(fileName, rawData,
				this->GetWidth(), this->GetHeight(),
				LodePNGColorType::LCT_GREY);
		}
		else if (this->pf == ColorSpace::PixelFormat::RGB)
		{
			lodepng::encode(fileName, rawData,
				this->GetWidth(), this->GetHeight(),
				LodePNGColorType::LCT_RGB);
		}
		else if (this->pf == ColorSpace::PixelFormat::RGBA)
		{
			lodepng::encode(fileName, rawData,
				this->GetWidth(), this->GetHeight(),
				LodePNGColorType::LCT_RGBA);
		}


	}

	if constexpr (std::is_same<T, uint8_t>::value == false)
	{
		delete rawData;
		rawData = nullptr;
	}
}

//=================================================================================================
// Getters
//=================================================================================================


template <typename T>
ColorSpace::PixelFormat Image2d<T>::GetPixelFormat() const noexcept
{
	return this->pf;
}

template <typename T>
size_t Image2d<T>::GetChannelsCount() const noexcept
{
	return this->channelsCount;
}


template <typename T>
int Image2d<T>::GetWidth() const noexcept
{
	return this->dim.w;
}

template <typename T>
int Image2d<T>::GetHeight() const noexcept
{
	return this->dim.h;
}

template <typename T>
size_t Image2d<T>::GetPixelsCount() const noexcept
{
	return size_t(this->dim.w) * size_t(this->dim.h);
}

template <typename T>
const ImageDimension & Image2d<T>::GetDimension() const noexcept
{
	return this->dim;
}

/// <summary>
/// Convert 1D index to 2D [x, y]
/// </summary>
/// <param name="index"></param>
/// <param name="x"></param>
/// <param name="y"></param>
template <typename T>
void Image2d<T>::GetPositionFromIndex(size_t index, int & x, int & y) const noexcept
{
	x = static_cast<int>(index % this->GetWidth());
	y = static_cast<int>(index / this->GetWidth());
}

/// <summary>
/// Convert 2D index [x, y] to 1D index
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
template <typename T>
size_t Image2d<T>::GetIndexFromPosition(int x, int y) const noexcept
{
	return size_t(x) + size_t(y) * size_t(this->dim.w);
}

/// <summary>
/// Get indices of 8-ring neighborhood of center pixel
/// Layout:
/// 0 | 1 | 2
/// 3 | x | 4
/// 5 | 6 | 7
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
template <typename T>
std::array<size_t, 8> Image2d<T>::GetNeighborsIndicesFromPosition(int x, int y) const noexcept
{
	int x1m = (x < 1) ? 0 : x - 1;
	int x1p = (x + 1 >= this->GetWidth()) ? this->GetWidth() - 1 : x + 1;

	int y1m = (y < 1) ? 0 : y - 1;
	int y1p = (y + 1 >= this->GetHeight()) ? this->GetHeight() - 1 : y + 1;

	std::array<size_t, 8> res;

	res[0] = this->GetIndexFromPosition(x1m, y1m);
	res[1] = this->GetIndexFromPosition(x1m, y);
	res[2] = this->GetIndexFromPosition(x1m, y1p);

	res[3] = this->GetIndexFromPosition(x, y1m);
	res[4] = this->GetIndexFromPosition(x, y1p);

	res[5] = this->GetIndexFromPosition(x1p, y1m);
	res[6] = this->GetIndexFromPosition(x1p, y);
	res[7] = this->GetIndexFromPosition(x1p, y1p);

	return res;
}

/// <summary>
/// Find min / max values for a given channelIndex
/// min / max are outputs
/// </summary>
/// <param name="min"></param>
/// <param name="max"></param>
/// <param name="channelIndex"></param>
template <typename T>
void Image2d<T>::FindMinMax(size_t channelIndex, T & min, T & max) const
{
#ifdef ENABLE_SIMD
	if constexpr (std::is_same<T, float>::value)
	{
		if (this->GetChannelsCount() == 1)
		{
			this->FindMinMaxSimd(min, max);
			return;
		}
	}
#endif

	const T * tmp = this->GetPixelStart(0);
	min = tmp[channelIndex];
	max = tmp[channelIndex];
	
	size_t len = this->GetPixelsCount();
	for (size_t i = 1; i < len; i++)
	{
		const T * tmp = this->GetPixelStart(i);

		max = std::max(tmp[channelIndex], max);
		min = std::min(tmp[channelIndex], min);		
	}
}

/// <summary>
/// Get pointer to the first element of pixel
/// </summary>
/// <param name="index"></param>
/// <returns></returns>
template <typename T>
const T * Image2d<T>::GetPixelStart(size_t index) const
{
	return &this->data[index * this->GetChannelsCount()];
}

template <typename T>
T * Image2d<T>::GetPixelStart(size_t index)
{
	return &this->data[index * this->GetChannelsCount()];
}


/// <summary>
/// Get pixel at position [x, y]
/// Position can be outside of image bounds and based on border mode,
/// it is updates to be in image bounds
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="border"></param>
/// <returns></returns>
template <typename T>
const T * Image2d<T>::GetPixelStartWithBorder(int x, int y, ImageUtils::BorderMode border) const
{	
	int h1 = this->dim.h - 1;
	int w1 = this->dim.w - 1;

	if (border == ImageUtils::BorderMode::WRAP)
	{
		y = (y < 0) ? h1 : (y > h1) ? 0 : y;
		x = (x < 0) ? w1 : (x > w1) ? 0 : x;
	}
	else if (border == ImageUtils::BorderMode::ENLARGE)
	{
		//reflect
		//[-2] [-1] [0] [1] [2]
		//put [-1] on [0]
		//put [-2] on [1]

		y = (y < 0) ? -y - 1 : (y > h1) ? (2 * this->dim.h - y - 1) : y;
		x = (x < 0) ? -x - 1 : (x > w1) ? (2 * this->dim.w - x - 1) : x;
	}
	else
	{
		//clamp
		y = (y < 0) ? 0 : (y > h1) ? h1 : y;
		x = (x < 0) ? 0 : (x > w1) ? w1 : x;
	}

	return this->GetPixelStart(x, y);
}

/// <summary>
/// Get pointer to the first byte of pixel
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param> 
/// <returns></returns>
template <typename T>
const T * Image2d<T>::GetPixelStart(int x, int y) const
{
	return this->GetPixelStart(this->GetIndexFromPosition(x, y));
}

template <typename T>
T * Image2d<T>::GetPixelStart(int x, int y)
{
	return this->GetPixelStart(this->GetIndexFromPosition(x, y));
}

template <typename T>
const T * Image2d<T>::operator[](size_t index) const
{
	return this->GetPixelStart(index);
}

template <typename T>
T * Image2d<T>::operator[](size_t index)
{
	return this->GetPixelStart(index);
}

template <typename T>
const std::vector<T> & Image2d<T>::GetData() const noexcept
{
	return this->data;
}

template <typename T>
std::vector<T> & Image2d<T>::GetData() noexcept
{
	return this->data;
}

//=================================================================================================
// Methods
//=================================================================================================

template <typename T>
void Image2d<T>::RunGauss(int radius, double weight, double strength)
{
	ImageGausFilter<T> gauss(radius, weight, strength);
	gauss.SetInput(*this);
	gauss.SetOutput(*this);
	gauss.RunSeparable();
	//gauss.RunBruteForce();
}

template <typename T>
void Image2d<T>::ForEachPixel(std::function<void(T *, size_t)> callback)
{
	size_t len = this->GetPixelsCount();

	for (size_t i = 0; i < len; i++)
	{
		T * tmp = this->GetPixelStart(i);
		callback(tmp, i);	
	}
	
}

template <typename T>
void Image2d<T>::ForEachPixel(std::function<void(const T *, size_t)> callback) const
{
	size_t len = this->GetPixelsCount();

	for (size_t i = 0; i < len; i++)
	{
		const T * tmp = this->GetPixelStart(i);
		callback(tmp, i);	
	}	
}

template <typename T>
void Image2d<T>::ForEachPixelInNeighborhood(size_t index, const NeighborhoodKernel & k,
	std::function<void(T *, int, int)> callback)
{
	int x, y;
	this->GetPositionFromIndex(index, x, y);
	this->ForEachPixelInNeighborhood(x, y, k, callback);
}


template <typename T>
void Image2d<T>::ForEachPixelInNeighborhood(size_t index, const NeighborhoodKernel & k,
	std::function<void(const T *, int, int)> callback) const
{
	int x, y;
	this->GetPositionFromIndex(index, x, y);
	this->ForEachPixelInNeighborhood(x, y, k, callback);
}

template <typename T>
void Image2d<T>::ForEachPixelInNeighborhood(int x, int y, const NeighborhoodKernel & k,
	std::function<void(T *, int, int)> callback)
{
	int yy_start, yy_stop;
	int xx_start, xx_stop;

	ImageUtils::CalcKernelStartStop(x, this->GetWidth() - 1, k, xx_start, xx_stop);
	ImageUtils::CalcKernelStartStop(y, this->GetHeight() - 1, k, yy_start, yy_stop);

	for (int yy = yy_start; yy <= yy_stop; yy++)
	{
		for (int xx = xx_start; xx <= xx_stop; xx++)
		{
			T * tmp = this->GetPixelStart(xx, yy);
			callback(tmp, xx, yy);
		}
	}
}

template <typename T>
void Image2d<T>::ForEachPixelInNeighborhood(int x, int y, const NeighborhoodKernel & k,
	std::function<void(const T *, int, int)> callback) const
{
	int yy_start, yy_stop;
	int xx_start, xx_stop;

	ImageUtils::CalcKernelStartStop(x, this->GetWidth() - 1, k, xx_start, xx_stop);
	ImageUtils::CalcKernelStartStop(y, this->GetHeight() - 1, k, yy_start, yy_stop);

	for (int yy = yy_start; yy <= yy_stop; yy++)
	{
		for (int xx = xx_start; xx <= xx_stop; xx++)
		{
			const T * tmp = this->GetPixelStart(xx, yy);
			callback(tmp, xx, yy);
		}
	}
}


template <typename T>
void Image2d<T>::ForEachPixelInLayeredNeighborhood(size_t index, const NeighborhoodKernel & k,
	std::function<void(T *, int, int, int)> callback)
{
	int x, y;
	this->GetPositionFromIndex(index, x, y);
	this->ForEachPixelInLayeredNeighborhood(x, y, k, callback);
}

/// <summary>
/// Iterate image layer by layer based on neighborhood size starting at pixel [x, y]
/// For each ¨processed pixel, calback is called (px, x, y, layer)
/// Layer is iterated clock-wise starting at top-left corner
/// (Note: If layer pixels would overflow the image, layer is not fully iterated,
/// but only pixels within the image are processed)
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="k"></param>
/// <param name="callback"></param>
template <typename T>
void Image2d<T>::ForEachPixelInLayeredNeighborhood(int x, int y, const NeighborhoodKernel & k,
	std::function<void(T *, int, int, int)> callback)
{	
	int layer = 0;
	int xx, yy;

	while (layer <= k.radius)
	{
		yy = y - layer;
		if ((yy > 0) && (yy < this->GetHeight()))
		{
			int start = std::max(x - layer, 0);
			int end = std::min(x + layer, this->GetWidth() - 1);

			for (xx = start; xx <= end; xx++)
			{
				T * tmp = this->GetPixelStart(xx, yy);
				callback(tmp, xx, yy, layer);
			}
		}

		xx = x + layer;
		if ((xx > 0) && (xx < this->GetWidth()))
		{
			int start = std::max(y - layer + 1, 0);
			int end = std::min(y + layer - 1, this->GetHeight() - 1);

			for (yy = start; yy <= end; yy++)
			{
				T * tmp = this->GetPixelStart(xx, yy);
				callback(tmp, xx, yy, layer);				
			}
		}

		yy = y + layer;
		if ((yy > 0) && (yy < this->GetHeight()))
		{
			int start = std::max(x - layer, 0);
			int end = std::min(x + layer, this->GetWidth() - 1);

			for (xx = end; xx >= start; xx--)
			{
				T * tmp = this->GetPixelStart(xx, yy);
				callback(tmp, xx, yy, layer);
			}
		}

		xx = x - layer;
		if ((xx > 0) && (xx < this->GetWidth()))
		{
			int start = std::max(y - layer + 1, 0);
			int end = std::min(y + layer - 1, this->GetHeight() - 1);

			for (yy = end; yy >= start; yy--)
			{
				T * tmp = this->GetPixelStart(xx, yy);
				callback(tmp, xx, yy, layer);
			}
		}

		layer++;
	}
}

/// <summary>
/// Calculaet absolute value of image
/// Only supported for float data type, since uint8_t is always positive
/// Use SIMD if enabled
/// </summary>
template <typename T>
void Image2d<T>::Abs()
{
	if constexpr (std::is_same<T, uint8_t>::value)
	{
		//not supported for uint8_t since it is already positive
		return;
	}
	else
	{
		size_t len8 = 0;

#ifdef ENABLE_SIMD
		len8 = this->data.size() - this->data.size() % MM256_ELEMENT_COUNT;

		for (size_t i = 0; i < len8; i += MM256_ELEMENT_COUNT)
		{
			auto v = _mm256_loadu_ps(this->data.data() + i);
			v = _my_mm256_abs_ps(v);
			_mm256_storeu_ps(this->data.data() + i, v);
		}

#endif

		for (size_t i = len8; i < this->data.size(); i++)
		{
			this->data[i] = std::abs(this->data[i]);
		}
	}
}

/// <summary>
/// Combine current image and input img
/// for each pixel, callback is called with user defined combination
/// Result can be written to a in callback (a = this)
/// </summary>
/// <param name="img"></param>
/// <param name="callback"></param>
template <typename T>
void Image2d<T>::Combine(const Image2d<T> & img, std::function<void(T *, const T *, size_t size)> callback)
{
	size_t len = this->GetPixelsCount();
	for (size_t i = 0; i < len; i++)
	{
		T * a = this->GetPixelStart(i);
		const T * b = img.GetPixelStart(i);
		callback(a, b, this->GetChannelsCount());
	}
}

/// <summary>
/// Map data to interval [start, end]
/// </summary>
/// <param name="start"></param>
/// <param name="end"></param>
template <typename T>
void Image2d<T>::MapToInterval(T start, T end)
{
#ifdef ENABLE_SIMD
	if constexpr (std::is_same<T, float>::value)
	{
		if (this->GetChannelsCount() == 1)
		{
			this->MapToIntervalSimd(start, end);
			return;
		}
	}
#endif

	size_t len = this->GetPixelsCount();
	
	T * maxVals = new T[this->GetChannelsCount()];
	T * minVals = new T[this->GetChannelsCount()];
			
	for (size_t c = 0; c < this->GetChannelsCount(); c++)
	{
		this->FindMinMax(c, minVals[c], maxVals[c]);
		
		//both values are the same - map to middle of interval
		if (minVals[c] == maxVals[c])
		{
			minVals[c] = 0;
			maxVals[c] = 1;
		}
	}
			
	for (size_t i = 0; i < len; i++)
	{
		T * tmp = this->GetPixelStart(i);
		
		for (size_t c = 0; c < this->GetChannelsCount(); c++)
		{
			float mappedVal = MathUtils::MapRange<float>(minVals[c], maxVals[c], start, end, tmp[c]);
			
			tmp[c] = static_cast<T>(mappedVal);			
		}
	}

	delete[] maxVals;
	delete[] minVals;
	maxVals = nullptr;
	minVals = nullptr;
	
}

#ifdef ENABLE_SIMD
template <>
void Image2d<uint8_t>::FindMinMaxSimd(uint8_t & min, uint8_t & max) const
{
	//not currently supported
}

template <>
void Image2d<uint8_t>::MapToIntervalSimd(uint8_t start, uint8_t end)
{
    //not currently supported
}

template<>
void Image2d<float>::FindMinMaxSimd(float & min, float & max) const
{
    size_t len8 = this->data.size() - this->data.size() % MM256_ELEMENT_COUNT;


    __m256 maxS = _mm256_loadu_ps(this->data.data());
    __m256 minS = maxS;

    for (size_t i = 8; i < len8; i += MM256_ELEMENT_COUNT)
    {
        auto v = _mm256_loadu_ps(this->data.data() + i);
        maxS = _mm256_max_ps(maxS, v);
        minS = _mm256_min_ps(minS, v);
    }

    max = _my_mm256_simple_hmax(maxS);
    min = _my_mm256_simple_hmin(minS);

    for (size_t i = len8; i < this->data.size(); i++)
    {
        max = std::max(max, this->data[i]);
        min = std::min(min, this->data[i]);
    }

    if (max == min)
    {
        min = 0;
        max = 1;
    }
}

/// <summary>
/// Map data to interval [start, end] using SIMD vectorization
/// </summary>
/// <param name="start"></param>
/// <param name="end"></param>
template <>
void Image2d<float>::MapToIntervalSimd(float start, float end)
{
	size_t len8 = this->data.size() - this->data.size() % MM256_ELEMENT_COUNT;

	//============================================================
	// Find max & min
	float min, max;
	this->FindMinMaxSimd(min, max);
	
	//============================================================

	float d = (end - start) / (max - min);

	auto minS = _mm256_set1_ps(min);	
	auto startS = _mm256_set1_ps(start);
	
	auto dS = _mm256_set1_ps(d);

	for (size_t i = 0; i < len8; i += MM256_ELEMENT_COUNT)
	{
		//MapRange
		//start + (s - min) * d

		auto s = _mm256_loadu_ps(this->data.data() + i);

		s = _mm256_sub_ps(s, minS);
		s = _mm256_mul_ps(s, dS);
		s = _mm256_add_ps(s, startS);

		_mm256_storeu_ps(this->data.data() + i, s);
	}

	for (size_t i = len8; i < this->data.size(); i++)
	{
		this->data[i] = start + (this->data[i] - min) * d;
	}

}


#endif

/// <summary>
/// Append image to the right of the current image
/// Appended image must have same number of channels as the current one
/// </summary>
/// <param name="img"></param>
template <typename T>
void Image2d<T>::AppendRight(const Image2d<T> & img)
{
	if (this->channelsCount != img.channelsCount)
	{
		MY_LOG_ERROR("Number of channels of append image is not same");
		return;
	}

	ImageDimension newDim;
	newDim.w = dim.w + img.GetWidth();
	newDim.h = std::max(dim.h, img.GetHeight());
		
	size_t len = size_t(newDim.w) * size_t(newDim.h) * this->channelsCount;

	std::vector<T> channelData;
	channelData.resize(len);

	//size_t newIndex = 0;

	//copy old image to the new array
	for (size_t y = 0; y < size_t(dim.h); y++)
	{
		std::copy(this->data.data() + ((0 + y * dim.w) * this->channelsCount),
			this->data.data() + ((dim.w + y * dim.w) * this->channelsCount),
			channelData.data() + ((0 + y * newDim.w) * this->channelsCount));

		/*
		for (size_t x = 0; x < dim.w; x++)
		{
			size_t oldIndex = (x + y * dim.w) * img.channelsCount;
			size_t newIndex = (x + y * newDim.w) * this->channelsCount;			
			for (size_t c = 0; c < this->channelsCount; c++)
			{
				channelData[newIndex + c] = this->data[oldIndex + c];
			}
		}
		*/
	}

	//copy appended image to the new array
	for (size_t y = 0; y < size_t(img.dim.h); y++)
	{
		std::copy(img.data.data() + ((0 + y * img.dim.w) * this->channelsCount),
			img.data.data() + ((img.dim.w + y * img.dim.w) * this->channelsCount),
			channelData.data() + ((dim.w + y * newDim.w) * this->channelsCount));
		
		/*
		for (size_t x = 0; x < img.dim.w; x++)
		{
			size_t oldIndex = (x + y * img.dim.w) * img.channelsCount;
			size_t newIndex = ((dim.w + x) + (y) * newDim.w) * this->channelsCount;			
			for (size_t c = 0; c < img.channelsCount; c++)
			{
				channelData[newIndex + c] = img.data[oldIndex + c];
			}
		}
		*/
	}
	
	this->dim = newDim;
	this->data = std::move(channelData);
}


/// <summary>
/// Helper method to convert gaus filtered data
/// to image data based on image type and channels count
/// </summary>
/// <param name="tmp"></param>
/// <param name="index"></param>
/// <returns></returns>
template <typename T>
void Image2d<T>::SetValue(double tmp, size_t channel, size_t index)
{
	T * val = this->GetPixelStart(index);
	val[channel] = ImageUtils::clamp_cast<T>(tmp);
}

template <typename T>
void Image2d<T>::SetValue(double tmp, size_t channel, int x, int y)
{
	T * val = this->GetPixelStart(x, y);
	val[channel] = ImageUtils::clamp_cast<T>(tmp);
}

template <typename T>
void Image2d<T>::SetValue(T value, size_t channel, size_t index)
{
	T * val = this->GetPixelStart(index);
	val[channel] = value;
}

template <typename T>
void Image2d<T>::SetValue(T value, size_t channel, int x, int y)
{
	T * val = this->GetPixelStart(x, y);
	val[channel] = value;
}

template <typename T>
void Image2d<T>::Clear(T clearValue)
{	
	std::fill(data.begin(), data.end(), clearValue);
}

template <typename T>
void Image2d<T>::SwapChannels(size_t c0, size_t c1)
{
	if ((c0 >= this->channelsCount) || (c1 >= this->channelsCount))
	{
		return;
	}

	size_t len = this->GetPixelsCount();
	
	for (size_t i = 0; i < len; i++)
	{
		T * val = this->GetPixelStart(i);

		T tmp = val[c0];
		val[c0] = val[c1];
		val[c1] = tmp;		
	}
}

template <typename T>
void Image2d<T>::AddChannels(size_t count)
{
	size_t len = this->GetPixelsCount();

	std::vector<T> d;
	d.resize(len * (this->channelsCount + count), 0);

	for (size_t i = 0; i < len; i++)
	{		
		const T * val = this->GetPixelStart(i);

		size_t newIndex = i * (this->channelsCount + count);
		for (size_t c = 0; c < this->channelsCount; c++)
		{
			d[newIndex + c] = val[c];
		}
	}

	this->channelsCount += count;
	this->pf = ColorSpace::PixelFormat::NONE;
}


//==================================================================================
// OpenCV only related section
//==================================================================================

#ifdef HAVE_OPENCV

/// <summary>
/// Get OpenCV Mat representation.
/// Note: 
/// It only holds pointer to the data, so Image2d instance must exist
/// while working with OpenCV Mat
/// </summary>
/// <returns></returns>
template <typename T>
cv::Mat Image2d<T>::CreateOpenCVLightCopy()
{
	int cvFormat = 0;

	if constexpr (std::is_same<T, uint8_t>::value)
	{
		if (this->pf == ColorSpace::PixelFormat::GRAY) cvFormat = CV_8UC1;
		else if (this->pf == ColorSpace::PixelFormat::RGB) cvFormat = CV_8UC3;
		else if (this->pf == ColorSpace::PixelFormat::RGBA) cvFormat = CV_8UC4;
		else
		{
			MY_LOG_ERROR("Incorrect pixel format for uint8_t image");
		}
	}
	else
	{
		//RGB and RGBA image can be float with values in range 0 - 1

		if (this->pf == ColorSpace::PixelFormat::GRAY) cvFormat = CV_32FC1;
		else if (this->pf == ColorSpace::PixelFormat::RGB) cvFormat = CV_32FC3;
		else if (this->pf == ColorSpace::PixelFormat::RGBA) cvFormat = CV_32FC4;
		else if (this->pf == ColorSpace::PixelFormat::XYZ) cvFormat = CV_32FC3;
		else if (this->pf == ColorSpace::PixelFormat::CIE_LUV) cvFormat = CV_32FC3;
		else if (this->pf == ColorSpace::PixelFormat::HSV) cvFormat = CV_32FC3;
		else
		{
			MY_LOG_ERROR("Incorrect pixel format for float image");
		}
	}

	return cv::Mat(static_cast<int>(this->GetHeight()),
		static_cast<int>(this->GetWidth()),
		cvFormat,
		this->data.data());
}

template <typename T>
cv::Mat Image2d<T>::CreateOpenCVDeepCopy()
{
	int cvFormat = 0;

	if constexpr (std::is_same<T, uint8_t>::value)
	{
		if (this->pf == ColorSpace::PixelFormat::GRAY) cvFormat = CV_8UC1;
		else if (this->pf == ColorSpace::PixelFormat::RGB) cvFormat = CV_8UC3;
		else if (this->pf == ColorSpace::PixelFormat::RGBA) cvFormat = CV_8UC4;
		else
		{
			MY_LOG_ERROR("Incorrect pixel format for uint8_t image");
		}
	}
	else
	{
		//RGB and RGBA image can be float with values in range 0 - 1

		if (this->pf == ColorSpace::PixelFormat::GRAY) cvFormat = CV_32FC1;
		else if (this->pf == ColorSpace::PixelFormat::RGB) cvFormat = CV_32FC3;
		else if (this->pf == ColorSpace::PixelFormat::RGBA) cvFormat = CV_32FC4;
		else if (this->pf == ColorSpace::PixelFormat::XYZ) cvFormat = CV_32FC3;
		else if (this->pf == ColorSpace::PixelFormat::CIE_LUV) cvFormat = CV_32FC3;
		else if (this->pf == ColorSpace::PixelFormat::HSV) cvFormat = CV_32FC3;
		else
		{
			MY_LOG_ERROR("Incorrect pixel format for float image");
		}
	}

	cv::Mat m(static_cast<int>(this->GetHeight()),
		static_cast<int>(this->GetWidth()),
		cvFormat);
	memcpy(m.data, this->data.data(), this->data.size() * sizeof(T));
	
	return m;
}

template <typename T>
Image2d<T> & Image2d<T>::operator=(const cv::Mat & cvMat)
{
	this->dim.w = cvMat.cols;
	this->dim.h = cvMat.rows;
	this->channelsCount = cvMat.channels();
	
	data.clear();

	if (cvMat.isContinuous())
	{
		data.assign((T*)cvMat.data, (T*)cvMat.data + cvMat.total() * channelsCount);
	}
	else
	{
		for (int i = 0; i < cvMat.rows; i++)
		{
			data.insert(data.end(), cvMat.ptr<T>(i), cvMat.ptr<T>(i) + cvMat.cols);
		}
	}

	return *this;
}

template <typename T>
void Image2d<T>::RunOpenCV(std::function<void(Image2d<T> * img, cv::Mat &)> cvCallback)
{
	cv::Mat cvMat = this->CreateOpenCVLightCopy();
	cvCallback(this, cvMat);
}

template <typename T>
void Image2d<T>::RunWithResultOpenCV(std::function<cv::Mat(Image2d<T> * img, cv::Mat &)> cvCallback)
{
	cv::Mat cvMat = this->CreateOpenCVLightCopy();
	auto res = cvCallback(this, cvMat);
	*this = res;
}

template <typename T>
void Image2d<T>::ShowImageOpenCV(bool waitForKey)
{

	cv::Mat cvMat = this->CreateOpenCVLightCopy();

	const char* cv_window = "OpenCV image";

	cv::namedWindow(cv_window, cv::WINDOW_NORMAL);
	cv::resizeWindow(cv_window, 
		std::min(this->GetWidth(), 1680), 
		std::min(this->GetHeight(), 1050));
	cv::imshow(cv_window, cvMat);
	if (waitForKey)
	{
		cv::waitKey();
	}
}
#endif


//=================================================================================================

template Image2d<float> Image2d<float>::CreateAs() const;
template Image2d<float> Image2d<uint8_t>::CreateAs() const;
template Image2d<uint8_t> Image2d<float>::CreateAs() const;
template Image2d<uint8_t> Image2d<uint8_t>::CreateAs() const;

template Image2d<float> Image2d<float>::CreateAsMapped(float start, float end) const;
template Image2d<float> Image2d<uint8_t>::CreateAsMapped(float start, float end) const;
template Image2d<uint8_t> Image2d<float>::CreateAsMapped(uint8_t start, uint8_t end) const;
template Image2d<uint8_t> Image2d<uint8_t>::CreateAsMapped(uint8_t start, uint8_t end) const;

template Image2d<float> Image2d<uint8_t>::CreateEmpty() const;
template Image2d<uint8_t> Image2d<uint8_t>::CreateEmpty() const;
template Image2d<float> Image2d<float>::CreateEmpty() const;
template Image2d<uint8_t> Image2d<float>::CreateEmpty() const;

template class Image2d<uint8_t>;
template class Image2d<float>;
