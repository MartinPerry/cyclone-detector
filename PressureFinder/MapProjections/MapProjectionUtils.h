#ifndef MAP_PROJECTION_UTILS_H
#define MAP_PROJECTION_UTILS_H

#include <vector>
#include <array>


#include "./MapProjectionStructures.h"
#include "./ProjectionInfo.h"

namespace Projections
{

	struct ProjectionUtils
	{
		
		static void ComputeAABB(const std::vector<Coordinate>& c, 
			Coordinate& min, Coordinate& max);
			
		static void ComputeAABB(const Coordinate& center, MyRealType dist,
			Coordinate& min, Coordinate& max);

		static void ComputeAABB_LessAccurateNearPoles(const Coordinate& center, MyRealType dist,
			Coordinate& min, Coordinate& max);

		static Coordinate CalcEndPointShortest(const Coordinate & start, const AngleValue & bearing, MyRealType dist);
		static Coordinate CalcEndPointDirect(const Coordinate & start, const AngleValue & bearing, MyRealType dist);
		static double Distance(const Coordinate & from, const Coordinate & to);
		
		static double CalcArea(const std::vector<Coordinate> & pts);
		
		template <typename PixelType, typename Projection>
		static double CalcArea(const std::vector<Pixel<PixelType>> & pxs, const Projection * from)
		{
			if (pxs.size() <= 2)
			{
				return 0.0;
			}

			std::vector<Coordinate> pts;
			for (auto & px : pxs)
			{
				auto gps = from->template  ProjectInverse<PixelType, false>(px);
				pts.push_back(gps);
			}

			return CalcArea(pts);
		}

		static std::array<Latitude, 2> EarthLatitudeRange(Latitude lat, double earthRadius, double distance);
		static std::array<Longitude, 2> EarthLongitudeRange(Latitude lat, Longitude lng, double earthRadius, int distance);
		static double CalcEarthRadiusAtLat(Latitude latitude);

        inline static MyRealType cot(MyRealType x) { return 1.0 / std::tan(x); };
        inline static MyRealType sec(MyRealType x) { return 1.0 / std::cos(x); };
        inline static MyRealType sinc(MyRealType x) { return std::sin(x) / x; };
        inline static MyRealType sgn(MyRealType x) { return (x < 0) ? -1 : (x > 0); };
          


	};

};

#endif
