#include "MapProjectionStructures.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <limits>
#include <algorithm>
#include <string.h>
#include <errno.h>

#ifndef MY_LOG_ERROR
#	define MY_LOG_ERROR(...) printf(__VA_ARGS__)
#endif

#ifdef HAVE_NEON
#	include "./Simd/neon_math_float.h"
#endif

using namespace Projections;

const double ProjectionConstants::PI = double(std::acos(-1));
const double ProjectionConstants::PI_4 = double(0.25) * ProjectionConstants::PI;
const double ProjectionConstants::PI_2 = double(0.5) * ProjectionConstants::PI;
const double ProjectionConstants::E = double(std::exp(1.0));
const double ProjectionConstants::EARTH_RADIUS = double(6371);


/// <summary>
/// Create Coordinate from cartexian [x, y, z]. 
/// Input is assumed to be in left-handed coordinate system
/// IS HAVE_NEON is defined - calculate 4 coordinates at once
/// 
/// Source:
/// https://vvvv.org/blog/polar-spherical-and-geographic-coordinates
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="z"></param>
/// <returns></returns>
std::array<Coordinate, 4> Coordinate::CreateFromCartesianLHSystem(
	const std::array<double, 4> & vx,
	const std::array<double, 4> & vy,
	const std::array<double, 4> & vz)
{
	std::array<Coordinate, 4> res;

#ifdef HAVE_NEON
	float32x4_t x =
	{
		static_cast<float32_t>(vx[0]),
		static_cast<float32_t>(vx[1]),
		static_cast<float32_t>(vx[2]),
		static_cast<float32_t>(vx[3])
};

	float32x4_t y =
	{
		static_cast<float32_t>(vy[0]),
		static_cast<float32_t>(vy[1]),
		static_cast<float32_t>(vy[2]),
		static_cast<float32_t>(vy[3])
	};

	float32x4_t z =
	{
		static_cast<float32_t>(vz[0]),
		static_cast<float32_t>(vz[1]),
		static_cast<float32_t>(vz[2]),
		static_cast<float32_t>(vz[3])
	};


	float32x4_t sum = vmulq_f32(x, x);
	sum = vaddq_f32(sum, vmulq_f32(y, y));
	sum = vaddq_f32(sum, vmulq_f32(z, z));

#if defined(__aarch64__) || defined(__arm64__) || (defined(vdivq_f32) && defined(vsqrtq_f32))
	float32x4_t radius = vsqrtq_f32(sum);
	float32x4_t rInv = vdivq_f32(vmovq_n_f32(1), radius);
#else
	float32x4_t radius = vmulq_f32(vrsqrteq_f32(sum), sum);
	float32x4_t rInv = vrecpeq_f32(radius); //estimate 1.0 / radius
#endif	

	float32x4_t lat = my_asin_f32(vmulq_f32(y, rInv));
	float32x4_t lon = my_atan2_f32(x, vmulq_n_f32(z, -1.0f));

	float resLat[4];
	float resLon[4];
	vst1q_f32(resLat, lat);
	vst1q_f32(resLon, lon);
	
	res[0] = Coordinate(Longitude::rad(resLon[0]), Latitude::rad(resLat[0]));
	res[1] = Coordinate(Longitude::rad(resLon[1]), Latitude::rad(resLat[1]));
	res[2] = Coordinate(Longitude::rad(resLon[2]), Latitude::rad(resLat[2]));
	res[3] = Coordinate(Longitude::rad(resLon[3]), Latitude::rad(resLat[3]));
#else
	res[0] = Coordinate::CreateFromCartesianLHSystem(vx[0], vy[0], vz[0]);
	res[1] = Coordinate::CreateFromCartesianLHSystem(vx[1], vy[1], vz[1]);
	res[2] = Coordinate::CreateFromCartesianLHSystem(vx[2], vy[2], vz[2]);
	res[3] = Coordinate::CreateFromCartesianLHSystem(vx[3], vy[3], vz[3]);
#endif
	return res;
}

/// <summary>
/// Create Coordinate from cartexian [x, y, z]. 
/// Input is assumed to be in left-handed coordinate system
/// 
/// Source:
/// https://vvvv.org/blog/polar-spherical-and-geographic-coordinates
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="z"></param>
/// <returns></returns>
Coordinate Coordinate::CreateFromCartesianLHSystem(double x, double y, double z)
{
	double radius;
	return Coordinate::CreateFromCartesianLHSystem(x, y, z, &radius);
};

/// <summary>
/// Create Coordinate from cartexian [x, y, z]. 
/// Input is assumed to be in left-handed coordinate system
/// 
/// Return also radius of converted position
/// 
/// Source:
/// https://vvvv.org/blog/polar-spherical-and-geographic-coordinates
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="z"></param>
/// <param name="radius"></param>
/// <returns></returns>
Coordinate Coordinate::CreateFromCartesianLHSystem(double x, double y, double z, double * radius)
{
	*radius = std::sqrt(x * x + y * y + z * z);
	return Coordinate::CreateFromCartesianLHSystem(x, y, z, *radius);
}

/// <summary>
/// Create Coordinate from cartexian [x, y, z] and precomputed radius
/// Input is assumed to be in left-handed coordinate system
/// 
/// Source:
/// https://vvvv.org/blog/polar-spherical-and-geographic-coordinates
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="z"></param>
/// <param name="radius"></param>
/// <returns></returns>
Coordinate Coordinate::CreateFromCartesianLHSystem(double x, double y, double z, double radius)
{
	double lat = std::asin(y / radius);
	double lon = std::atan2(x, -z);
	return Coordinate(Longitude::rad(lon), Latitude::rad(lat));
}

/// <summary>
/// Convert cartesian 3D vector to (lat, lon) vector
/// in left-handed system
/// 
/// Based on:
/// https://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
/// (modified for left-handed system and for direct use of lat/lon)
/// </summary>
/// <param name="dx"></param>
/// <param name="dy"></param>
/// <param name="dz"></param>
/// <param name="precomp"></param>
/// <returns></returns>
Coordinate Coordinate::ConvertVectorFromCartesianLHSystem(double dx, double dy, double dz,
	const Coordinate::PrecomputedSinCos * precomp) const
{
	double rDif;
	return Coordinate::ConvertVectorFromCartesianLHSystem(dx, dy, dz, rDif, precomp);
}

/// <summary>
/// Convert cartesian 3D vector to (lat, lon) vector
/// in left-handed system
/// 
/// Return also rDif - difference vector of radius (= movement along the radius)
/// 
/// Based on:
/// https://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
/// (modified for left-handed system and for direct use of lat/lon)
/// </summary>
/// <param name="dx"></param>
/// <param name="dy"></param>
/// <param name="dz"></param>
/// <param name="rDif"></param>
/// <param name="precomp"></param>
/// <returns></returns>
Coordinate Coordinate::ConvertVectorFromCartesianLHSystem(double dx, double dy, double dz,
	double & rDif,
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

	rDif = cosLat * sinLon * dx + sinLat * dy - cosLat * cosLon * dz;
	double latDif = -sinLat * sinLon * dx + cosLat * dy + sinLat * cosLon * dz;
	double lonDif = cosLon * dx + sinLon * dz;

	return Coordinate(Longitude::rad(lonDif), Latitude::rad(latDif));
}

/// <summary>
/// Precompute sin/cos for lat/lon that can be used for repeated 
/// multiple computations
/// 
/// </summary>
/// <returns></returns>
Coordinate::PrecomputedSinCos Coordinate::PrecomputeSinCos() const
{
	Coordinate::PrecomputedSinCos tmp;

#ifdef HAVE_NEON
	//we can use neon only to calculate 2 sin/cos at once
	//other two possible calculations are not used
	float32x4_t v = { 
		static_cast<float32_t>(lat.rad()), 
		static_cast<float32_t>(lon.rad()),
		static_cast<float32_t>(0.0), 
		static_cast<float32_t>(0.0) 
	};

	float32x4_t ysin;
	float32x4_t ycos;
	my_sincos_f32(v, &ysin, &ycos);

	float values[4];
	vst1q_f32(values, ysin);
	tmp.sinLat = values[0];
	tmp.sinLon = values[1];

	vst1q_f32(values, ycos);
	tmp.cosLat = values[0];
	tmp.cosLon = values[1];

#else
	tmp.sinLat = std::sin(lat.rad());
	tmp.cosLat = std::cos(lat.rad());
	tmp.sinLon = std::sin(lon.rad());
	tmp.cosLon = std::cos(lon.rad());
#endif
	return tmp;
}

/// <summary>
/// Prepare precalculated sin/cos - faster for NEON because we can utilise all
/// calculations
/// </summary>
/// <param name="coords"></param>
/// <returns></returns>
std::array<Coordinate::PrecomputedSinCos, 4> Coordinate::PrecalcMultipleSinCos(
	const std::array<Coordinate, 4> & coords)
{
	std::array<Coordinate::PrecomputedSinCos, 4> res;

#ifdef HAVE_NEON
	for (size_t i = 0; i < 4; i += 2)
	{
		float32x4_t v = {
			static_cast<float32_t>(coords[i].lat.rad()),
			static_cast<float32_t>(coords[i].lon.rad()),
			static_cast<float32_t>(coords[i + 1].lat.rad()),
			static_cast<float32_t>(coords[i + 1].lon.rad())
		};

		float32x4_t ysin;
		float32x4_t ycos;
		my_sincos_f32(v, &ysin, &ycos);

		float values[4];
		vst1q_f32(values, ysin);
		res[i].sinLat = values[0];
		res[i].sinLon = values[1];
		res[i + 1].sinLat = values[2];
		res[i + 1].sinLon = values[3];

		vst1q_f32(values, ycos);
		res[i].cosLat = values[0];
		res[i].cosLon = values[1];
		res[i + 1].cosLat = values[2];
		res[i + 1].cosLon = values[3];
	}

#else
	for (size_t i = 0; i < 4; i++)
	{
		res[i].sinLat = std::sin(coords[i].lat.rad());
		res[i].cosLat = std::cos(coords[i].lat.rad());
		res[i].sinLon = std::sin(coords[i].lon.rad());
		res[i].cosLon = std::cos(coords[i].lon.rad());
	}
#endif

	return res;
}
