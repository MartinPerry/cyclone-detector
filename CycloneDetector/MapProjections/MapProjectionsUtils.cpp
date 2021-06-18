#include "./MapProjectionUtils.h"

#include "./MapProjectionStructures.h"

using namespace Projections;


void ProjectionUtils::ComputeAABB(const std::vector<Coordinate>& c,
	Coordinate& min, Coordinate& max)
{
	if (c.size() == 0)
	{
		return;
	}

	min = c[0];
	max = c[0];
	for (size_t i = 1; i < c.size(); i++)
	{
		if (c[i].lat.rad() < min.lat.rad()) min.lat = c[i].lat;
		if (c[i].lon.rad() < min.lon.rad()) min.lon = c[i].lon;


		if (c[i].lat.rad() > max.lat.rad()) max.lat = c[i].lat;
		if (c[i].lon.rad() > max.lon.rad()) max.lon = c[i].lon;

	}
}


/// <summary>
/// https://stackoverflow.com/questions/1689096/calculating-bounding-box-a-certain-distance-away-from-a-lat-long-coordinate-in-j
/// http://janmatuschek.de/LatitudeLongitudeBoundingCoordinates
/// </summary>
/// <param name="center"></param>
/// <param name="dist"></param>
/// <param name="min"></param>
/// <param name="max"></param>
void ProjectionUtils::ComputeAABB_LessAccurateNearPoles(const Coordinate& center, MyRealType dist,
	Coordinate& min, Coordinate& max)
{

	static const double MIN_LAT = -ProjectionConstants::PI_2;  // -PI/2
	static const double MAX_LAT = ProjectionConstants::PI_2;   //  PI/2
	static const double MIN_LON = -ProjectionConstants::PI; // -PI
	static const double MAX_LON = ProjectionConstants::PI;  //  PI

	// angular distance in radians on a great circle
	MyRealType radDist = dist / ProjectionConstants::EARTH_RADIUS;

	MyRealType minLat = center.lat.rad() - radDist;
	MyRealType maxLat = center.lat.rad() + radDist;

	double minLon, maxLon;
	if (minLat > MIN_LAT && maxLat < MAX_LAT) 
	{
		double deltaLon = asin(sin(radDist) / cos(center.lat.rad()));

		minLon = center.lon.rad() - deltaLon;
		if (minLon < MIN_LON) minLon += 2.0 * ProjectionConstants::PI;

		maxLon = center.lon.rad() + deltaLon;
		if (maxLon > MAX_LON) maxLon -= 2.0 * ProjectionConstants::PI;
	}
	else 
	{
		// a pole is within the distance
		minLat = std::max(minLat, MIN_LAT);
		maxLat = std::min(maxLat, MAX_LAT);
		minLon = MIN_LON;
		maxLon = MAX_LON;
	}

	min = Coordinate(Longitude::rad(minLon), Latitude::rad(minLat));
	max = Coordinate(Longitude::rad(maxLon), Latitude::rad(maxLat));
}

/// <summary>
/// Calculate end point based on shortest path (on real earth surface)
/// see: http://www.movable-type.co.uk/scripts/latlong.html
/// </summary>
/// <param name="start"></param>
/// <param name="bearing"></param>
/// <param name="dist"></param>
/// <returns></returns>
Coordinate ProjectionUtils::CalcEndPointShortest(const Coordinate & start,
	const AngleValue & bearing, MyRealType dist)
{

	MyRealType dr = dist / ProjectionConstants::EARTH_RADIUS;

	MyRealType sinLat = std::sin(start.lat.rad());
	MyRealType cosLat = std::cos(start.lat.rad());
	MyRealType sinDr = std::sin(dr);
	MyRealType cosDr = std::cos(dr);
	MyRealType sinBear = std::sin(bearing.rad());
	MyRealType cosBear = std::cos(bearing.rad());

	MyRealType sinEndLat = sinLat * cosDr + cosLat * sinDr * cosBear;
	MyRealType endLat = std::asin(sinEndLat);
	MyRealType y = sinBear * sinDr * cosLat;
	MyRealType x = cosDr - sinLat * sinLat;
	MyRealType endLon = start.lon.rad() + std::atan2(y, x);


	Coordinate end;
	end.lat = Latitude::rad(endLat);
	end.lon = Longitude::rad(endLon);
	end.lon.Normalize();

	return end;

}

/// <summary>
/// "Rhumb lines"
/// Calculate end point based on direct path (straight line between two points in projected earth)
/// </summary>
/// <param name="start"></param>
/// <param name="bearing"></param>
/// <param name="dist"></param>
/// <returns></returns>
Coordinate ProjectionUtils::CalcEndPointDirect(
	const Coordinate & start, const AngleValue & bearing, MyRealType dist)
{
	MyRealType dr = dist / ProjectionConstants::EARTH_RADIUS;

	MyRealType difDr = dr * std::cos(bearing.rad());
	MyRealType endLat = start.lat.rad() + difDr;

	// check for some daft bugger going past the pole, normalise latitude if so
	if (std::abs(endLat) > ProjectionConstants::PI_2) endLat = endLat > 0 ? ProjectionConstants::PI - endLat : -ProjectionConstants::PI - endLat;


	MyRealType projLatDif = std::log(std::tan(endLat / 2 + ProjectionConstants::PI_4) / std::tan(start.lat.rad() / 2 + ProjectionConstants::PI_4));
	MyRealType q = std::abs(projLatDif) > 10e-12 ? difDr / projLatDif : std::cos(start.lat.rad()); // E-W course becomes ill-conditioned with 0/0

	MyRealType difDrQ = dr * std::sin(bearing.rad()) / q;
	MyRealType endLon = start.lon.rad() + difDrQ;



	Coordinate end;
	end.lat = Latitude::rad(endLat);
	end.lon = Longitude::rad(endLon);
	end.lon.Normalize();

	return end;
}

/// <summary>
/// Calculate Haversine distance in km between from - to
/// 
/// /http://stackoverflow.com/questions/365826/calculate-distance-between-2-gps-coordinates
/// </summary>		
/// <param name="from"></param>
/// <param name="to"></param>
/// <returns></returns>
double ProjectionUtils::Distance(const Coordinate & from, const Coordinate & to)
{
	MyRealType dlong = to.lon.rad() - from.lon.rad();
	MyRealType dlat = to.lat.rad() - from.lat.rad();

	MyRealType a = std::pow(std::sin(dlat / 2.0), 2) + std::cos(from.lat.rad()) * std::cos(to.lat.rad()) * std::pow(std::sin(dlong / 2.0), 2);
	MyRealType c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
	MyRealType d = 6367 * c;

	if (dlong >= 3.14159265358979323846)
	{
		//we are going over 0deg meridian
		//distance maybe wrapped around the word - shortest path

		//split computation to [-lon, 0] & [0, lon]
		//which basically mean, subtract equator length => 40075km

		d = 40075.0 - d;
	}

	return d;
};

/// <summary>
/// Calculate area of polygon in meters
/// https://stackoverflow.com/questions/2861272/polygon-area-calculation-using-latitude-and-longitude-generated-from-cartesian-s
/// </summary>
/// <param name="pts"></param>
/// <returns></returns>
double ProjectionUtils::CalcArea(const std::vector<Coordinate> & pts)
{
	double area = 0.0;

	if (pts.size() > 2)
	{
		for (size_t i = 0; i < pts.size() - 1; i++)
		{
			Coordinate p1 = pts[i];
			Coordinate p2 = pts[i + 1];
			
			double lonDif = p2.lon.rad() - p1.lon.rad();
			area += lonDif * (2 + std::sin(p1.lat.rad()) + std::sin(p2.lat.rad()));
		}

		//last and first point to enclose polygon
		Coordinate p1 = pts[pts.size() - 1];
		Coordinate p2 = pts[0];

		double lonDif = p2.lon.rad() - p1.lon.rad();
		area += lonDif * (2 + std::sin(p1.lat.rad()) + std::sin(p2.lat.rad()));
		//-------


		//6378137 - earth radius in m
		area = area * 6378137 * 6378137 / 2;
	}

	return abs(area);
}



/// <summary>
/// Calculate latitude range based on earths radius at a given point
/// </summary>
/// <param name="lat"></param>
/// <param name="radius"></param>
/// <param name="distance"></param>
/// <returns></returns>
std::array<Latitude, 2> ProjectionUtils::EarthLatitudeRange(Latitude lat, double earthRadius, double distance) 
{
	// Estimate the min and max latitudes within distance of a given location.

	double angle = distance / earthRadius;
	double minlat = lat.rad() - angle;
	double maxlat = lat.rad() + angle;
	double rightangle = ProjectionConstants::PI / 2;
	// Wrapped around the south pole.
	if (minlat < -rightangle) 
	{
		double overshoot = -minlat - rightangle;
		minlat = -rightangle + overshoot;
		if (minlat > maxlat) {
			maxlat = minlat;
		}
		minlat = -rightangle;
	}
	// Wrapped around the north pole.
	if (maxlat > rightangle) 
	{
		double overshoot = maxlat - rightangle;
		maxlat = rightangle - overshoot;
		if (maxlat < minlat) {
			minlat = maxlat;
		}
		maxlat = rightangle;
	}
		
	return { Latitude::rad(minlat), Latitude::rad(maxlat) };
}

/// <summary>
/// Calculate longitude range based on earths radius at a given point
/// </summary>
/// <param name="lat"></param>
/// <param name="lng"></param>
/// <param name="earth_radius"></param>
/// <param name="distance"></param>
/// <returns></returns>
std::array<Longitude, 2> ProjectionUtils::EarthLongitudeRange(Latitude lat, Longitude lng, double earthRadius, int distance)
{
	// Estimate the min and max longitudes within distance of a given location.
	double radius = earthRadius * cos(lat.rad());

	double angle;
	if (radius > 0) 
	{
		angle = abs(distance / radius);
		angle = std::min(angle, ProjectionConstants::PI);
	}
	else 
	{
		angle = ProjectionConstants::PI;
	}
	double minlong = lng.rad() - angle;
	double maxlong = lng.rad() + angle;
	if (minlong < -ProjectionConstants::PI) 
	{
		minlong = minlong + ProjectionConstants::PI * 2;
	}
	if (maxlong > ProjectionConstants::PI) 
	{
		maxlong = maxlong - ProjectionConstants::PI * 2;
	}
	
	return { Longitude::rad(minlong), Longitude::rad(maxlong) };
}

/// <summary>
/// Calculate earth radius at given latitude
/// </summary>
/// <param name="latitude"></param>
/// <returns></returns>
double ProjectionUtils::CalcEarthRadiusAtLat(Latitude latitude)
{
	// Estimate the Earth's radius at a given latitude.
	// Default to an approximate average radius for the United States.
	double lat = latitude.rad();

	double x = cos(lat) / 6378137.0;
	double y = sin(lat) / (6378137.0 * (1 - (1 / 298.257223563)));

	//Make sure earth's radius is in km , not meters
	return (1 / (sqrt(x * x + y * y))) / 1000;
}


void ProjectionUtils::ComputeAABB(const Coordinate& center, MyRealType dist,
	Coordinate& min, Coordinate& max)
{	
	double radius = ProjectionUtils::CalcEarthRadiusAtLat(center.lat);
	auto retLats = ProjectionUtils::EarthLatitudeRange(center.lat, radius, dist);
	auto retLngs = ProjectionUtils::EarthLongitudeRange(center.lat, center.lon, radius, dist);

	min = Coordinate(retLngs[0], retLats[0]);
	max = Coordinate(retLngs[1], retLats[1]);
}
