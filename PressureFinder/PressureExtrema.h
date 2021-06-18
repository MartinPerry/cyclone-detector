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

//============================================================================

struct AreaThreshold
{
	enum class Unit { Pixels, Km2 };

	double value;
	Unit unit;
	Projections::IProjectionInfo* proj;

	constexpr static AreaThreshold CreateKm2(double value, Projections::IProjectionInfo* proj)
	{
		return AreaThreshold(value, AreaThreshold::Unit::Km2, proj);
	}

	constexpr static AreaThreshold CreatePixels(double value)
	{
		return AreaThreshold(value, AreaThreshold::Unit::Pixels, nullptr);
	}

	constexpr AreaThreshold() :
		value(0),
		unit(Unit::Pixels),
		proj(nullptr)
	{
	}

	constexpr AreaThreshold(double value, Unit unit, Projections::IProjectionInfo* proj) :
		value(value),
		unit(unit),
		proj(proj)
	{
	}

	constexpr AreaThreshold(const AreaThreshold & area) :
		value(area.value),
		unit(area.unit),
		proj(area.proj)
	{
	}

	AreaThreshold& operator=(AreaThreshold other)
	{		
		value = other.value;
		unit = other.unit;
		proj = other.proj;
		return *this;
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
	

	PressureExtrema(const Contour & c);

	bool IsValid() const;

	void AddPixel(int x, int y);
	
	void GetCenterIndex(int & cx, int & cy) const;

	int GetWidth() const;
	int GetHeight() const;

	double GetSizePixels() const;
	double GetSizeKm2(Projections::IProjectionInfo* proj) const;

	void ClearPixels();
	void Clear();


	void CreateAreas(PressureType type, const AreaThreshold& areaThresholdSize);
private:


	struct PressureArea
	{
		int minX = std::numeric_limits<int>::max();
		int minY = std::numeric_limits<int>::max();
		int maxX = 0;
		int maxY = 0;
		std::vector<ImageUtils::Pixel> pixels;

		void Clear()
		{
			minX = std::numeric_limits<int>::max();
			minY = std::numeric_limits<int>::max();
			maxX = 0;
			maxY = 0;
			pixels.clear();
		}

		void GetCenterIndex(int & cx, int & cy) const
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

	std::list<PressureArea> CreateAreasInner();
};



#endif
