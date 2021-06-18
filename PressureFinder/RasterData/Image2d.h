#ifndef IMAGE_2D_H
#define IMAGE_2D_H

struct NeighborhoodKernel;

#include <vector>
#include <functional>
#include <optional>

#if __has_include (<opencv2/imgproc.hpp>)
#	include <opencv2/imgproc.hpp>
#	define HAVE_OPENCV 1
#endif


#include "./ImageUtils.h"
#include "./ColorSpace.h"

template <typename T>
class Image2d 
{
public:
	
	static Image2d<T> CreateFromRawFile(int w, int h, ColorSpace::PixelFormat pf, const char * fileName);
	
	Image2d();		
	Image2d(const char * fileName);	
	Image2d(int w, int h, ColorSpace::PixelFormat pf);
	Image2d(int w, int h, const std::vector<T> & data, ColorSpace::PixelFormat pf);
	Image2d(int w, int h, const T * rawData, ColorSpace::PixelFormat pf);
	Image2d(int w, int h, const T ** rawData, ColorSpace::PixelFormat pf);
	Image2d(int w, int h, std::vector<T> && data, ColorSpace::PixelFormat pf);
	Image2d(const Image2d<T> & other);
	Image2d(Image2d<T> && other) noexcept;
#ifdef HAVE_OPENCV
	Image2d(const cv::Mat & cvMat);
#endif
	~Image2d();
	
	void Release() noexcept;

	Image2d<T> & operator=(const Image2d<T> & other);
	Image2d<T> & operator=(Image2d<T> && other) noexcept;

	template <typename V = T>
	Image2d<V> CreateEmpty() const;
	Image2d<T> CreateDeepCopy() const;
	Image2d<T> CreateFromChannel(size_t channelIndex) const;
	Image2d<T> CreateSubImage(int x, int y, const ImageDimension & size) const;
	
	template <typename V>
	Image2d<V> CreateAs() const;

	template <typename V>
	Image2d<V> CreateAsMapped(V start, V end) const;
	
#ifdef HAVE_OPENCV
	cv::Mat CreateOpenCVLightCopy();
	cv::Mat CreateOpenCVDeepCopy();
#endif
	
	void SetPixelFormat(ColorSpace::PixelFormat pf) noexcept;
	
	void Save(const char * fileName) const;

	ColorSpace::PixelFormat GetPixelFormat() const noexcept;
	size_t GetChannelsCount() const noexcept;
	

	int GetWidth() const noexcept;
	int GetHeight() const noexcept;
	size_t GetPixelsCount() const noexcept;
	const ImageDimension & GetDimension() const noexcept;
	
	void GetPositionFromIndex(size_t index, int & x, int & y) const noexcept;
	size_t GetIndexFromPosition(int x, int y) const noexcept;
	std::array<size_t, 8> GetNeighborsIndicesFromPosition(int x, int y) const noexcept;
	
	void FindMinMax(size_t channelIndex, T & min, T & max) const;

	const T * GetPixelStart(size_t index) const;
	const T * GetPixelStart(int x, int y) const;
	const T * GetPixelStartWithBorder(int x, int y, ImageUtils::BorderMode border) const;
	T * GetPixelStart(size_t index);
	T * GetPixelStart(int x, int y);
	const T * operator[](size_t index) const;
	T * operator[](size_t index);

	const std::vector<T> & GetData() const noexcept;
	std::vector<T> & GetData() noexcept;

	void RunGauss(int radius, double weight = 1.0, double strength = 1.0);
	
	void ForEachPixel(std::function<void(T *, size_t)> callback);
	void ForEachPixel(std::function<void(const T *, size_t)> callback) const;

	void ForEachPixelInNeighborhood(int x, int y, const NeighborhoodKernel & k,
		std::function<void(T *, int, int)> callback);
	void ForEachPixelInNeighborhood(int x, int y, const NeighborhoodKernel & k,
		std::function<void(const T *, int, int)> callback) const;
	void ForEachPixelInNeighborhood(size_t index, const NeighborhoodKernel & k,
		std::function<void(T *, int, int)> callback);
	void ForEachPixelInNeighborhood(size_t index, const NeighborhoodKernel & k,
		std::function<void(const T *, int, int)> callback) const;

	void ForEachPixelInLayeredNeighborhood(size_t index, const NeighborhoodKernel & k,
		std::function<void(T *, int, int, int)> callback);
	void ForEachPixelInLayeredNeighborhood(int x, int y, const NeighborhoodKernel & k,
		std::function<void(T *, int, int, int)> callback);
	


	void MapToInterval(T start, T end);	
	void Abs();	
	void Combine(const Image2d<T> & img, std::function<void(T *, const T *, size_t size)> callback);
	void AppendRight(const Image2d<T> & img);
	void SetValue(double tmp, size_t channel, size_t index);
	void SetValue(double tmp, size_t channel, int x, int y);
	void SetValue(T value, size_t channel, size_t index);
	void SetValue(T value, size_t channel, int x, int y);

	void Clear(T clearValue = T(0));
	void SwapChannels(size_t c0, size_t c1);
	void AddChannels(size_t count);
	
#ifdef HAVE_OPENCV
	Image2d<T> & operator=(const cv::Mat & cvMat);	
	void RunOpenCV(std::function<void(Image2d<T> * img, cv::Mat &)> cvCallback);
	void RunWithResultOpenCV(std::function<cv::Mat(Image2d<T> * img, cv::Mat &)> cvCallback);
	void ShowImageOpenCV(bool waitForKey = true);
#endif

    //https://stackoverflow.com/questions/3292795/how-to-declare-a-templated-struct-class-as-a-friend
    //We have to create Image2d friend with itself
    //so we can access protected members of float from uint8_t and vice versa
    template<typename> friend class Image2d;
    
protected:	
	ImageDimension dim;
	std::vector<T> data;
	ColorSpace::PixelFormat pf;	
	size_t channelsCount;

#ifdef ENABLE_SIMD
	void FindMinMaxSimd(T & min, T & max) const;
	void MapToIntervalSimd(T start, T end);
#endif
};


#endif
