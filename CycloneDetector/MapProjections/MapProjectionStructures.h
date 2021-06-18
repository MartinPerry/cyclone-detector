#ifndef MAP_PROJECTION_STRUCTURES_H
#define MAP_PROJECTION_STRUCTURES_H


#include <vector>
#include <string>
#include <tuple>
#include <array>

#include "GeoCoordinate.h"

#ifdef _MSC_VER
#	ifndef my_fopen 
#		define my_fopen(a, b, c) fopen_s(a, b, c)	
#	endif
#else
#	ifndef my_fopen 
#		define my_fopen(a, b, c) (*a = fopen(b, c))
#	endif
#endif

namespace Projections
{
	//================================================================================================
	//================================================================================================
	//================================================================================================


	enum class PROJECTION
    {
		UNKNOWN = -1,
		MERCATOR = 0,
		LAMBERT_CONIC = 1,
		EQUIRECTANGULAR = 2,
		WEB_MERCATOR = 3,
		POLAR_STEREOGRAPHICS = 4,
        MILLER = 5,
		GEOS = 6
	};


	enum STEP_TYPE
    {
		PIXEL_BORDER = 0,
		PIXEL_CENTER = 1
	};

	//================================================================================================
	//================================================================================================
	//================================================================================================


	template <typename PixelType = int,
    typename = typename std::enable_if<std::is_arithmetic<PixelType>::value, PixelType>::type>
    struct Pixel
    {
        PixelType x;
        PixelType y;
    };

	//================================================================================================
	//================================================================================================
	//================================================================================================


	struct Coordinate 
	{

		struct PrecomputedSinCos
		{
			double sinLat;
			double cosLat;
			double sinLon;
			double cosLon;
		};

		Longitude lon;
		Latitude lat;

		Coordinate() {};
		Coordinate(Longitude lon, Latitude lat) :
			lon(lon), lat(lat) {};

		PrecomputedSinCos PrecomputeSinCos() const;

		static std::array<Coordinate::PrecomputedSinCos, 4> PrecalcMultipleSinCos(
			const std::array<Coordinate, 4> & coords);

		static std::array<Coordinate, 4> CreateFromCartesianLHSystem(
			const std::array<double, 4> & x,
			const std::array<double, 4> & y,
			const std::array<double, 4> & z);

		static Coordinate CreateFromCartesianLHSystem(double x, double y, double z);
		static Coordinate CreateFromCartesianLHSystem(double x, double y, double z, double * radius);
		static Coordinate CreateFromCartesianLHSystem(double x, double y, double z, double radius);

		template <typename T = std::tuple<double, double, double>>
		T ConvertToCartesianLHSystem(double radius,
			const Coordinate::PrecomputedSinCos * precomp = nullptr) const;

		Coordinate ConvertVectorFromCartesianLHSystem(double dx, double dy, double dz,
			const Coordinate::PrecomputedSinCos * precomp = nullptr) const;

		Coordinate ConvertVectorFromCartesianLHSystem(double dx, double dy, double dz, double & radius,
			const Coordinate::PrecomputedSinCos * precomp = nullptr) const;

		template <typename T = std::tuple<double, double, double>>
		T ConvertVectorToCartesianLHSystem(const Longitude & lonDif, const Latitude & latDif,
			const Coordinate::PrecomputedSinCos * precomp = nullptr) const;
	};

		
	
	//================================================================================================
	//================================================================================================
	//================================================================================================


	struct ProjectionFrame
	{
		Coordinate min;
		Coordinate max;
	   
		MyRealType w; //current frame width
		MyRealType h; //current frame height

		MyRealType wPadding;
		MyRealType hPadding;
		MyRealType wAR; //width AR
		MyRealType hAR; //height AR

		MyRealType projPrecomX;
		MyRealType projPrecomY;

		STEP_TYPE stepType;


		MyRealType repeatNegCount;
		MyRealType repeatPosCount;

		MyRealType ww; //width of single full projected world
		MyRealType hh; //height of single full projected world

		ProjectionFrame() : 
			min(Coordinate(Longitude::deg(-180.0), Latitude::deg(-90.0))),
			max(Coordinate(Longitude::deg(180.0), Latitude::deg(90.0))),
			w(0),
			h(0),
			wPadding(0),
			hPadding(0),
			wAR(1),
			hAR(1),
			projPrecomX(0),
			projPrecomY(0),
			stepType(STEP_TYPE::PIXEL_CENTER),
			repeatNegCount(0),
			repeatPosCount(0),
			ww(0),
			hh(0)
		{}
	};

	//================================================================================================
	//================================================================================================
	//================================================================================================


    struct ProjectionConstants
    {
        static const MyRealType PI;
        static const MyRealType PI_4;
        static const MyRealType PI_2;
        static const MyRealType E;
        static const MyRealType EARTH_RADIUS;
    };
    
	//================================================================================================
	//================================================================================================
	//================================================================================================

	/// <summary>
	/// Convert latitude / longitude to left-handed cartesian coordinate system
	/// Radius is used as scaled earth radius
	/// 
	/// Source:
	/// https://vvvv.org/blog/polar-spherical-and-geographic-coordinates
	/// </summary>
	/// <param name="radius"></param>
	/// <param name="precomp"></param>
	/// <returns></returns>
	template <typename T>
	T Coordinate::ConvertToCartesianLHSystem(double radius, const Coordinate::PrecomputedSinCos * precomp) const
	{
		double sinLat, cosLat, sinLon, cosLon;

		if (precomp != nullptr)
		{
			sinLat = precomp->sinLat;
			cosLat = precomp->cosLat;
			sinLon = precomp->sinLon;
			cosLon = precomp->cosLon;
		}
		else
		{
			Coordinate::PrecomputedSinCos tmp = this->PrecomputeSinCos();
			sinLat = tmp.sinLat;
			cosLat = tmp.cosLat;
			sinLon = tmp.sinLon;
			cosLon = tmp.cosLon;
		}

		double x = radius * cosLat * sinLon;
		double y = radius * sinLat;
		double z = -radius * cosLat * cosLon;

		return { x, y, z };
	}

	/// <summary>
	/// Convert (lat, lon) vector to cartesian coordinate vector
	/// in left-handed system
	/// 
	/// Based on:
	/// https://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
	/// (modified for left-handed system and for direct use of lat/lon)
	/// </summary>
	/// <param name="lonDif"></param>
	/// <param name="latDif"></param>
	/// <param name="precomp"></param>
	/// <returns></returns>
	template <typename T>
	T Coordinate::ConvertVectorToCartesianLHSystem(const Longitude & lonDif, const Latitude & latDif,
		const Coordinate::PrecomputedSinCos * precomp) const
	{
		double sinLat, cosLat, sinLon, cosLon;

		if (precomp != nullptr)
		{
			sinLat = precomp->sinLat;
			cosLat = precomp->cosLat;
			sinLon = precomp->sinLon;
			cosLon = precomp->cosLon;
		}
		else
		{
			Coordinate::PrecomputedSinCos tmp = this->PrecomputeSinCos();
			sinLat = tmp.sinLat;
			cosLat = tmp.cosLat;
			sinLon = tmp.sinLon;
			cosLon = tmp.cosLon;
		}

		//full calculation if vector contains also movement in radius direction
		//double dx = cosLat * sinLon  * (r) - sinLat * sinLon * (latDif.rad()) + cosLon * (lonDif.rad());
		//double dy = sinLat * (r) + cosLat * (latDif.rad());
		//double dz = -cosLat * cosLon * (r) + sinLat * cosLon  * (latDif.rad()) + sinLon * (lonDif.rad());

		double dx = -sinLat * sinLon * (latDif.rad()) + cosLon * (lonDif.rad());
		double dy = cosLat * (latDif.rad());
		double dz = sinLat * cosLon  * (latDif.rad()) + sinLon * (lonDif.rad());

		return { dx, dy, dz };
	}
};

#endif
