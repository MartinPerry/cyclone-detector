#ifndef EQUIRECTANGULAR_H
#define EQUIRECTANGULAR_H

#include <cmath>

#include "../GeoCoordinate.h"
#include "../ProjectionInfo.h"
#include "../MapProjectionStructures.h"


namespace Projections
{

	/// <summary>
	/// Based on:
	/// https://en.wikipedia.org/wiki/Equirectangular_projection
	/// </summary>
	class Equirectangular : public ProjectionInfo<Equirectangular>
	{
	public:

		static const bool INDEPENDENT_LAT_LON = true; //can Lat / Lon be computed separatly. To compute one, we dont need the other
		static const bool ORTHOGONAL_LAT_LON = true; //lat / lon is orthogonal to each otjer

		Equirectangular() : Equirectangular(Longitude::deg(0.0)) {}
		Equirectangular(const Longitude & lonCentralMeridian) :
			ProjectionInfo(PROJECTION::EQUIRECTANGULAR),
			lonCentralMeridian(lonCentralMeridian),
			standardParallel(Latitude::deg(0.0)),
			cosStandardParallel(std::cos(standardParallel.rad()))
		{ }

		friend class ProjectionInfo<Equirectangular>;

	protected:

		const Longitude lonCentralMeridian;
		const Latitude standardParallel;
		const MyRealType cosStandardParallel;

		const char* GetNameInternal() const
		{
			return "Equirectangular";
		}

		ProjectedValue ProjectInternal(const Coordinate & c) const
		{
			return {
				(c.lon.rad() - lonCentralMeridian.rad()) * cosStandardParallel,
				c.lat.rad() - standardParallel.rad()
			};
		};

		ProjectedValueInverse ProjectInverseInternal(MyRealType x, MyRealType y) const
		{			
			return {
				Latitude::rad(y / cosStandardParallel + lonCentralMeridian.rad()),
				Longitude::rad(x + standardParallel.rad())
			};
		};

	};
}

#endif
