#include "./SimplificationDouglasPeucker.h"

#include <algorithm>


#include "../../Math/Vector2.h"

#ifdef MYMATH_NAMESPACE
using namespace MYMATH_NAMESPACE;
#endif

#ifdef GRAPHICS_NAMESPACE
using namespace GRAPHICS_NAMESPACE;
#endif

using namespace d2;

SimplificationDouglasPeucker::SimplificationDouglasPeucker(const std::vector<Vector2d> & pts) :
	pts(pts)
{
}

double SimplificationDouglasPeucker::DistanceSquared(const Vector2d & a, const Vector2d & b) const
{
	double x = a.x - b.x;
	double y = a.y - b.y;
	return (x * x) + (y * y);
}

double SimplificationDouglasPeucker::PerpedincularDistanceSquared(const Vector2d & a, const Vector2d & b, const Vector2d & p) const
{
	// Return minimum distance between line segment ab and point p
	const double len2 = DistanceSquared(a, b);  // i.e. |b-a|^2 -  avoid a sqrt
	if (len2 == 0.0)
	{
		return DistanceSquared(p, a);   // a == b case
	}
	// Consider the line extending the segment, parameterized as a + t (b - a).
	// We find projection of point p onto the line. 
	// It falls where t = [(p-a) . (b-a)] / |b-a|^2
	// We clamp t from [0,1] to handle points outside the segment vw.
	
	// dot(p - a, b - a)
	const double dx1 = p.x - a.x;
	const double dy1 = p.y - a.y;

	const double dx2 = b.x - a.x;
	const double dy2 = b.y - a.y;

	const double tmp = (dx1 * dx2 + dy1 * dy2) / len2;
	//

	const double t = std::max(0.0, std::min(1.0, tmp));

	const double projectionX = a.x + t * dx2;  // Projection falls on the segment
	const double projectionY = a.y + t * dy2;  // Projection falls on the segment
	return DistanceSquared(p, { projectionX, projectionY });
}

/// <summary>
/// Simplify line with Douglas-Peucker
/// </summary>
/// <param name="threshold"></param>
/// <returns></returns>
PolyLine SimplificationDouglasPeucker::Simplify(double threshold) const
{	

	//opravit sebe-protinani

	PolyLine simpleLine;

	//there is not enough points to create simplification triangles
	//create manually simplified lines from the same input points
	if (this->pts.size() < 3)
	{
		for (auto & p : this->pts)
		{
			simpleLine.AddPoint(p);
		}

		return simpleLine;
	}

	std::vector<int> result;

	result.push_back(0); //startin point is always present	
	this->Simplify(0, static_cast<int>(pts.size()) - 1, threshold, result);

	if (!(pts[result.back()] == pts[static_cast<int>(pts.size()) - 1]))
	{
		result.push_back(static_cast<int>(pts.size()) - 1);
	}

	for (auto p : result)
	{
		simpleLine.AddPoint(pts[p]);
	}
	
	

	//printf("Before: %d\n", pts.size());
	//printf("After: %d\n", simpleLine.Size());

	return simpleLine;
}


void SimplificationDouglasPeucker::Simplify(int startIndex, int endIndex, double tresshold, 
	std::vector<int> & result) const
{
	if (endIndex - startIndex < 2)
	{
		//no points between - add end point
		//result.push_back(startIndex);
		result.push_back(endIndex);
		return;
	}

	//find max distance from straight line
	double max = 0;
	int maxIndex = 0;

	for (int i = startIndex + 1; i < endIndex; i++)
	{
		double d = PerpedincularDistanceSquared(pts[startIndex], pts[endIndex], pts[i]);

		if (d > max)
		{
			maxIndex = i;
			max = d;
		}
	}


	if (max > tresshold)
	{
		//Recursive call
		Simplify(startIndex, maxIndex, tresshold, result);
		Simplify(maxIndex, endIndex, tresshold, result);

	}
	else
	{
		//if (pts[result.back()] == pts[startIndex])
		{
			//if last index value in result == value at start index in polyline
			//insert only "end" point
			result.push_back(endIndex);
		}		
	}
}