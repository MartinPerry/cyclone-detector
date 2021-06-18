#ifndef IMAGE_UTILS_H
#define IMAGE_UTILS_H

template <typename T>
class Image2d;

#include <stdint.h>
#include <vector>
#include <list>
#include <array>
#include <functional>

//=========================================================================

/// <summary>
/// Image dimension structure
/// Holds width and height
/// </summary>
struct ImageDimension
{
	int w;
	int h;

	/// <summary>
	/// Create resized image with given scale
	/// </summary>
	/// <param name="scale"></param>
	/// <returns></returns>
	ImageDimension CreateResized(double scale) const
	{
		return { static_cast<int>(w * scale), static_cast<int>(h * scale) };
	};

	/// <summary>
	/// Create new dimension and keep aspect ratio
	/// new width is calculated from given new height and AR of current size
	/// </summary>
	/// <param name="newHeight"></param>
	/// <returns></returns>
	ImageDimension CreateResizedWidthAr(int newHeight) const
	{
		double ar = double(w) / h;
		return { static_cast<int>(newHeight * ar), newHeight };
	};

	/// <summary>
	/// Create new dimension and keep aspect ratio
	/// new height is calculated from given new width and AR of current size
	/// </summary>
	/// <param name="newWidth"></param>
	/// <returns></returns>
	ImageDimension CreateResizedHeightAr(int newWidth) const
	{
		double ar = double(h) / w;
		return { newWidth, static_cast<int>(newWidth * ar) };
	};
};

//=========================================================================

struct NeighborhoodKernel
{
	enum class Connectivity 
	{
		CON_4 = 4,
		CON_8 = 8
	};

	static NeighborhoodKernel CreateFromSize(int size);
	static NeighborhoodKernel CreateFromRadius(int r);

	int radius; //filter radius ("border around center pixel")

	int GetSize() const;
	int GetPixelArea() const;

private:
	NeighborhoodKernel(int radius);

};

//=========================================================================

class ImageUtils 
{
public:
	enum class BorderMode 
	{ 
		CLAMP = 0, 
		WRAP = 1, 
		ENLARGE = 2,
		ZERO = 3	//border is set to 0
	};

	struct CornerDetection
	{
		//for corner detection itself
		int gausRadius;
		double gausWeight;
		double gausStrength;

		float kappa; //valid only for CornerHarris function
		float(*cornerFunctionCallback)(float, float, float, float);

		//used if we look for best corners
		NeighborhoodKernel k;
		float percentage;	
		float threshold;

		CornerDetection() : 
			gausRadius(1),
			gausWeight(1.0),
			gausStrength(1.0),
			kappa(0.04f),
			cornerFunctionCallback(ImageUtils::CornerHarrisFunction),			
			k(NeighborhoodKernel::CreateFromSize(0)),
			percentage(0),
			threshold(0)
		{
		}
	};

	struct Pixel 
	{
		int x;
		int y;

		Pixel() : x(0), y(0) {}
		Pixel(int x, int y) : x(x), y(y) {}
	};

	template <typename T>
	struct TileInfo 
	{
		int x;
		int y;
		Image2d<T> img;
	};

	static std::array<uint8_t, 3> CreateRandomColor();

	template <typename T>
	static void DrawLine(Image2d<T> & input, const T * value,
		int x0, int y0, int x1, int y1);
	
	static void ProcessLinePixels(int x0, int y0, int x1, int y1,
		std::function<void(int x, int y)> pixelCallback);

	template <typename T>
	static size_t FloodFill(Image2d<T> & input, int x, int y,
		T newVal, 
		std::function<bool(T, T, size_t)> sameAreaCallback,
		NeighborhoodKernel::Connectivity con);

	template <typename T>
	static size_t FloodFill(Image2d<T> & input, int x, int y,
		T newVal, T valDifThreshold, 
		NeighborhoodKernel::Connectivity con);

	template <typename T>
	static std::list<std::vector<Pixel>> Vectorize(const Image2d<T> & input, T bgValue, 
		double minLength2 = 0, 
		bool joinLines = false);

	template <typename T>
	static Image2d<float> SDF(const Image2d<T> & input, T bgValue);

	
	static void CalcKernelStartStop(double v0, int maxVal,
		const NeighborhoodKernel & k,
		int & start, int & end);
	
	static void CalcKernelStartStop(int v0, int maxVal,
		const NeighborhoodKernel & k,
		int & start, int & end);

	template <typename T>
	static void CalcDerivatives(const Image2d<T> & input, int x, int y,
		float * dx, float * dy, float * dxx, float * dyy, float * dxy);

	template <typename T>
	static Image2d<float> DetectCornersHarris(const Image2d<T> & input, const CornerDetection & info,
		std::vector<size_t> * bestCorners);

	template <typename T, typename V>
	static T clamp_cast(const V & val);

	static std::vector<size_t> FindMaxValues(const Image2d<float> & input, const NeighborhoodKernel & k, float percents);

	template <typename T>
	static std::vector<TileInfo<T>> CreateTiles(const Image2d<T> & input, int tileW, int tileH, bool allSameSize = true);


	static float CornerHarrisFunction(float a, float b, float c, float k);
	static float CornerShiTomasiFunction(float a, float b, float c, float k);

private:
	static const int INSIDE = 0; // 0000
	static const int LEFT = 1;   // 0001
	static const int RIGHT = 2;  // 0010
	static const int BOTTOM = 4; // 0100
	static const int TOP = 8;    // 1000

	static const float INF;// = 1E20;

	static int ComputeOutCode(int x, int y, int w, int h);
	static void JoinVectorizedLines(int w, int h, std::list<std::vector<Pixel>> & lines);

	static void AddVectorizedLine(std::vector<ImageUtils::Pixel> & newLine, 
		int sinceLastSplitCount,
		double lineLen,
		double minLength2,
		std::list<std::vector<Pixel>> & dest);

	static void SDFTransformToDistance(Image2d<float> & img);
	static void SDFCalcEdt(const float* f, float* d, float* z, uint16_t* w, int n);
};


/// <summary>
/// cast input val to output
/// if output is uint8_t, 
/// cast is with a clamp to interval [0,255]
/// </summary>
/// <param name="val"></param>
/// <returns></returns>
template <typename T, typename V>
T ImageUtils::clamp_cast(const V & val)
{
	if constexpr (std::is_same<T, uint8_t>::value)
	{
		return (double(val) > 255.0) ? static_cast<T>(255) :
			(double(val) < 0.0) ? static_cast<T>(0) : static_cast<T>(val);
	}
	else
	{
		return static_cast<T>(val);
	}
}

#endif
