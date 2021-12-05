#ifndef PRESSURE_EXTREMA_H
#define PRESSURE_EXTREMA_H


namespace Projections
{
	class IProjectionInfo;
}


#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include <list>
#include <vector>

#include "./Math/Vector2.h"

#include "./Graphics/2d/Polygon.h"
#include "./RasterData/ImageUtils.h"

#include "./MapProjections/MapProjectionStructures.h"

//============================================================================
struct PressureThreshold
{
public:
	enum class Unit { Pixels, Km2, Km };

	double value;
	Unit unit;
	Projections::IProjectionInfo* proj;


	constexpr PressureThreshold(const PressureThreshold& area) :
		value(area.value),
		unit(area.unit),
		proj(area.proj)
	{
	}

	PressureThreshold& operator=(PressureThreshold other)
	{
		value = other.value;
		unit = other.unit;
		proj = other.proj;
		return *this;
	}

protected:

	constexpr PressureThreshold(double value, Unit unit, Projections::IProjectionInfo* proj) :
		value(value),
		unit(unit),
		proj(proj)
	{
	}
};

struct AreaThreshold : public PressureThreshold
{
public:
	constexpr static AreaThreshold CreateKm2(double value, Projections::IProjectionInfo* proj)
	{
		return AreaThreshold(value, PressureThreshold::Unit::Km2, proj);
	}

	constexpr static AreaThreshold CreatePixels(double value)
	{
		return AreaThreshold(value, PressureThreshold::Unit::Pixels, nullptr);
	}

	constexpr AreaThreshold(double value, Unit unit, Projections::IProjectionInfo* proj) :
		PressureThreshold(value, unit, proj)
	{
	}
};

struct DistanceThreshold : public PressureThreshold
{
public:
	constexpr static DistanceThreshold CreateKm(double value, Projections::IProjectionInfo* proj)
	{
		return DistanceThreshold(value, PressureThreshold::Unit::Km, proj);
	}

	constexpr static DistanceThreshold CreatePixels(double value)
	{
		return DistanceThreshold(value, PressureThreshold::Unit::Pixels, nullptr);
	}

	constexpr DistanceThreshold(double value, Unit unit, Projections::IProjectionInfo* proj) :
		PressureThreshold(value, unit, proj)
	{
	}
};

//============================================================================

/// <summary>
/// Contour info 
/// </summary>
struct Contour
{
	double t; //contour value
	d2::Polygon c; //contour itself, points in pixels
	Contour * parent; //parent contour
	
	int childCount; //number of inner contours
	std::vector<Contour *> childs; //list of inner contours
	
	bool isExtrema; //is contour detected as "extrema" (it encloses extrema)

	Contour(double t, const d2::Polygon & c) :
		t(t),
		c(c),
		parent(nullptr),
		childCount(0),
		isExtrema(false)
	{}

	
	Contour(double t, d2::Polygon && c) :
		t(t),
		c(std::move(c)),
		parent(nullptr),
		childCount(0),
		isExtrema(false)
	{}	
};

//============================================================================

/// <summary>
/// Extrema info
/// </summary>
class PressureExtrema
{
public:	
	enum class PressureType 
	{
		NONE = -1,
		LOW = 0,
		HIGH = 1
	};

	const Contour & contour;
	
	PressureType type;
	MyMath::Vector2f center;
	double centerValue;
	

	PressureExtrema(const Contour & c);

	bool IsValid() const;

	void AddPixel(int x, int y);
	
	Projections::Coordinate GetCenterGps(Projections::IProjectionInfo* proj) const;
	void GetAabbCenterIndex(int& cx, int& cy) const;

	int GetWidth() const;
	int GetHeight() const;

	double GetSizePixels() const;
	double GetSizeKm2(Projections::IProjectionInfo* proj) const;

	void ClearPixels();
	void Clear();

	bool HasPixels() const;
	bool CheckArea(const AreaThreshold& areaThresholdSize) const;
	void InitCenter(PressureType type, const Image2d<float>& rawData);

private:


	struct PressureArea
	{
		int minX = std::numeric_limits<int>::max();
		int minY = std::numeric_limits<int>::max();
		int maxX = 0;
		int maxY = 0;
		std::vector<ImageUtils::Pixel> pixels;

		ImageUtils::Pixel minValuePx;
		ImageUtils::Pixel maxValuePx;
		double minValue = std::numeric_limits<double>::max();
		double maxValue = std::numeric_limits<double>::lowest();

		void Clear()
		{
			minValuePx.x = -1;
			minValuePx.y = -1;

			maxValuePx.x = -1;
			maxValuePx.y = -1;

			minValue = std::numeric_limits<double>::max();
			maxValue = std::numeric_limits<double>::lowest();

			minX = std::numeric_limits<int>::max();
			minY = std::numeric_limits<int>::max();
			maxX = 0;
			maxY = 0;
			pixels.clear();
		}

		/// <summary>
		/// Get center index as center of AABB - fast buf not accurate
		/// </summary>
		/// <param name="cx"></param>
		/// <param name="cy"></param>
		void GetAabbCenterIndex(int& cx, int& cy) const
		{
			cx = minX + (maxX - minX) / 2;
			cy = minY + (maxY - minY) / 2;
		}


		int GetSize() const
		{
			return ((maxX - minX) + 1) * ((maxY - minY) + 1);
		}
	};
	
	std::optional<PressureArea> bestPressure;

	std::vector<ImageUtils::Pixel> px;

	int minX = std::numeric_limits<int>::max();
	int minY = std::numeric_limits<int>::max();
	int maxX = 0;
	int maxY = 0;

	void GetPositionFromIndex(int index, int & x, int & y) const;

	int GetIndexFromPosition(int x, int y) const;

	bool TestAreaPixel(size_t ind, std::unordered_set<int> & pixels,
		std::unordered_map<int, PressureArea *> & result);

	std::list<PressureArea> CreateAreasInner(const Image2d<float>& rawData);
};



#endif
