#include "./Polygon.h"

#include <cmath>
#include <algorithm>

#include "../../Math/Vector2.h"

#ifdef MYMATH_NAMESPACE
using namespace MYMATH_NAMESPACE;
#endif


#ifdef GRAPHICS_NAMESPACE
using namespace GRAPHICS_NAMESPACE;
#endif

using namespace d2;

Polygon Polygon::CreateFromPolyline(const PolyLine & l)
{
	Polygon p;
	p.bb = l.bb;
	p.pts = l.pts;
	if (l.pts.size() > 1)
	{
		p.AddPoint(l.pts[0]);
	}
	return p;
}

Polygon::Polygon() noexcept
{
}

Polygon::Polygon(PolyLine && line)
{
	this->bb = std::move(line.bb);
	this->pts = std::move(line.pts);
	if (this->pts.size() > 1)
	{
		auto firstPt = this->pts[0];
		this->AddPoint(firstPt);
	}
}

Polygon::Polygon(const Polygon & p)
{
	this->bb = p.bb;
	this->pts = p.pts;
}

/// <summary>
/// move ctor
/// </summary>
/// <param name="p"></param>
Polygon::Polygon(Polygon && p) noexcept	
{
	this->bb = std::move(p.bb);
	this->pts = std::move(p.pts);
}

/// <summary>
/// Copy assignment
/// </summary>
/// <param name="p"></param>
/// <returns></returns>
Polygon& Polygon::operator=(const Polygon& p)
{
	// Self-assignment detection
	if (&p == this)
	{
		return *this;
	}

	// Release any resource we're holding


	// Copy the resource
	this->bb = p.bb;
	this->pts = p.pts;

	return *this;
}

/// <summary>
/// Calc area of simple polygon
/// </summary>
/// <returns></returns>
double Polygon::CalcArea() const
{
	double a = 0.0, b = 0.0;
	for (size_t i = 0; i < this->pts.size() - 1; i++)
	{
		a += this->pts[i].x * this->pts[i + 1].y;
		b += this->pts[i].y * this->pts[i + 1].x;
	}
	return std::abs((a - b) / 2.0);
}

/// <summary>
/// Test if any part of polygon is in box
/// </summary>
/// <param name="box"></param>
/// <returns></returns>
bool Polygon::IsInBox(const Aabb<Vector2d> & box) const
{
	bool inBox = PolyLine::IsInBox(box);
	if (inBox)
	{
		return true;
	}

	//test enclosing line between last and first point
	const Vector2d & curPt = this->pts.back();
	const Vector2d & nextPt = this->pts[0];

	int curCode = box.GetCode(curPt);
	int nextCode = box.GetCode(nextPt);

	if (!(curCode | nextCode))
	{
		// bitwise OR is 0: both points inside window
		//at least two points inside - accepth
		return true;
	}
	else if (curCode & nextCode)
	{
		// bitwise AND is not 0: both points are outside window; we know nothing
		return false;
	}

	//intersection
	return true;
}

/// <summary>
/// Test if point is inside complex polygon
/// http://alienryderflex.com/polygon/
/// </summary>
/// <param name="p"></param>
/// <returns></returns>
bool Polygon::IsPointInside(const Vector2d & p) const
{
	if (this->bb.IsInside(p) == false)
	{
		return false;
	}

	size_t j = this->pts.size() - 1;
	bool  oddNodes = false;

	for (size_t i = 0; i < this->pts.size(); i++)
	{
		double iX = this->pts[i].x;
		double iY = this->pts[i].y;

		double jX = this->pts[j].x;
		double jY = this->pts[j].y;

		if (((iY < p.y && jY >= p.y) || (jY < p.y && iY >= p.y)) && (iX <= p.x || jX <= p.x))
		{
			oddNodes ^= (iX + (p.y - iY) / (jY - iY)*(jX - iX) < p.x);
		}
		j = i;
	}

	/*
	for (size_t i = 0; i < this->pts.size(); i++)
	{
		if (this->pts[i].y < p.y && this->pts[j].y >= p.y || this->pts[j].y < p.y && this->pts[i].y >= p.y)
		{
			if (this->pts[i].x + (p.y - this->pts[i].y) / (this->pts[j].y - this->pts[i].y)*(this->pts[j].x - this->pts[i].x) < p.x)
			{
				oddNodes = !oddNodes;
			}
		}
		j = i;
	}
	*/

	return oddNodes;
}

/// <summary>
/// Calc centroid of polygon
/// https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
/// </summary>
/// <returns></returns>
Vector2d Polygon::CalcCentroid() const
{
	Vector2d centroid = Vector2d(0, 0);

	double signedArea = 0.0;
	double x0 = 0.0; // Current vertex X
	double y0 = 0.0; // Current vertex Y
	double x1 = 0.0; // Next vertex X
	double y1 = 0.0; // Next vertex Y
	double a = 0.0;  // Partial signed area

	// For all vertices except last
	size_t i = 0;
	for (i = 0; i < this->pts.size() - 1; i++)
	{
		x0 = this->pts[i].x;
		y0 = this->pts[i].y;
		x1 = this->pts[i + 1].x;
		y1 = this->pts[i + 1].y;
		a = x0 * y1 - x1 * y0;
		signedArea += a;
		centroid.x += (x0 + x1) * a;
		centroid.y += (y0 + y1) * a;
	}

	// Do last vertex
	x0 = this->pts[i].x;
	y0 = this->pts[i].y;
	x1 = this->pts[0].x;
	y1 = this->pts[0].y;
	a = x0 * y1 - x1 * y0;
	signedArea += a;
	centroid.x += (x0 + x1) * a;
	centroid.y += (y0 + y1) * a;

	signedArea *= 0.5;
	centroid.x /= (6 * signedArea);
	centroid.y /= (6 * signedArea);

	return centroid;
}

/// <summary>
/// Clamp polygon using combination of line clipping via Cohen-Sutherland
/// and Weiler–Atherton for polygon clipping
/// </summary>
/// <param name="head"></param>
/// <param name="box"></param>
/// <returns></returns>
std::list<Polygon> Polygon::Clamp(const Aabb<Vector2d> & box) const
{
	if (this->pts.size() == 0)
	{
		return std::list<Polygon>();
	}

	std::vector<ClampInfo> clampInfo;

	//add temporarily first point as last point to enclose polygon
	bool needAdd = false;
	if (this->pts[0] != this->pts.back())
	{
		needAdd = true;
		this->pts.push_back(this->pts.front());
	}

	std::list<PolyLine> clampedLines = this->ClampInternal(box, clampInfo);


	if (needAdd)
	{
		//remove temporary added point
		this->pts.pop_back();
	}


	if (clampedLines.size() == 0)
	{
		//data are empty - no clamping, but box can be fully inside polygon
		//we need to check this
		//=> test if corners are inside complex polygon
		//http://alienryderflex.com/polygon/



		auto p0 = box.GetCorner(0);
		bool oddNodes = this->IsPointInside(p0);

		if (oddNodes)
		{
			//one corner is inside and there was no clamping => all corners are isnide
			std::list<Polygon> polys;
			Polygon pol;
			pol.bb = box;
			pol.pts.push_back(p0);
			pol.pts.push_back(box.GetCorner(1));
			pol.pts.push_back(box.GetCorner(2));
			pol.pts.push_back(box.GetCorner(3));


			polys.push_back(std::move(pol));

			return polys;
		}

		//nothing in visible area
		return std::list<Polygon>();
	}

	if ((clampInfo[0].in == -1) && (clampInfo[0].out == -1))
	{
		std::list<Polygon> polys;
		for (PolyLine & l : clampedLines)
		{
			Polygon pol;
			pol.pts = std::move(l.pts);
			pol.bb = l.bb;
			polys.push_back(std::move(pol));
		}

		//no clamping - data are totaly in				
		return polys;
	}

	//========================================
	// Weiler–Atherton part of algorithm
	//========================================

	std::list<Polygon> result;
	Polygon pol;
	pol.bb = box;

	std::list<ClampPixel> clip = this->BuildClipWindow(clampedLines, clampInfo, box);

	std::list<PolyLine>::iterator pIt = clampedLines.begin();
	do
	{
		for (auto & px : pIt->pts)
		{
			pol.pts.push_back(px);
		}

		//find "exit" point in clip windows list
		auto wIt = std::find(clip.begin(), clip.end(), pol.pts.back());

		//circulation - if we have reached end => go to begin
		wIt = (std::next(wIt) == clip.end()) ? clip.begin() : std::next(wIt);

		//while next clip point is corner - add it
		while (wIt->corner)
		{
			pol.pts.push_back(*wIt);

			//move to next point
			//circulation - if we have reached end => go to begin
			wIt = (std::next(wIt) == clip.end()) ? clip.begin() : std::next(wIt);
		}


		clampedLines.erase(pIt);
		pIt = clampedLines.end();

		//find "entry" point in polygon list
		for (auto it = clampedLines.begin(); it != clampedLines.end(); it++)
		{
			if (it->pts.front() == *wIt)
			{
				pIt = it;
				break;
			}
		}

		if (pIt == clampedLines.end())
		{
			//entry point not found
			//polygon is closed
			result.push_back(std::move(pol));
			pol.pts.clear();

			pIt = clampedLines.begin();
		}

	} while (!clampedLines.empty());


	return result;
}

/// <summary>
/// Build clipping window (AABB) list for Weiler-Atherton
/// </summary>
/// <param name="clampedLines"></param>
/// <param name="clampInfo"></param>
/// <param name="box"></param>
/// <returns></returns>
std::list<Polygon::ClampPixel> Polygon::BuildClipWindow(const std::list<PolyLine> & clampedLines,
	const std::vector<ClampInfo> & clampInfo, const Aabb<Vector2d> & box) const
{

	//build clipping window list
	std::list<ClampPixel> clipEdges[4];
	for (int i = 0; i < 4; i++)
	{
		//add starting corners
		clipEdges[i].emplace_back(box.GetCorner(i), true);
	}

	size_t i = 0;
	for (auto & line : clampedLines)
	{
		//add entry and exit points for each line
		//we have to push_front, because pushed point can be same as corner point
		//and in that case, sorting would result in pushed point being "last"
		//but we need "corner" as last after sorting
		clipEdges[clampInfo[i].in].emplace_front(line.pts.front(), false);
		clipEdges[clampInfo[i].out].emplace_front(line.pts.back(), false);
		i++;
	}

	//sort lines in CW order
	//corner is always "last"
	clipEdges[Aabb<Vector2d>::CLAMP_TOP].sort([](const Vector2d & a, const Vector2d & b) {
		return a.x > b.x;
		});
	clipEdges[Aabb<Vector2d>::CLAMP_LEFT].sort([](const Vector2d & a, const Vector2d & b) {
		return a.y > b.y;
		});
	clipEdges[Aabb<Vector2d>::CLAMP_BOTTOM].sort([](const Vector2d & a, const Vector2d & b) {
		return a.x < b.x;
		});
	clipEdges[Aabb<Vector2d>::CLAMP_RIGHT].sort([](const Vector2d & a, const Vector2d & b) {
		return a.y < b.y;
		});

	//join edges together
	std::list<ClampPixel> clip;
	clip.splice(clip.end(), clipEdges[Aabb<Vector2d>::CLAMP_TOP]);
	clip.splice(clip.end(), clipEdges[Aabb<Vector2d>::CLAMP_LEFT]);
	clip.splice(clip.end(), clipEdges[Aabb<Vector2d>::CLAMP_BOTTOM]);
	clip.splice(clip.end(), clipEdges[Aabb<Vector2d>::CLAMP_RIGHT]);

	return clip;
}

/// <summary>
/// Simplify Polygon as line with Visvalingam  algorithm
/// https://bost.ocks.org/mike/simplify/
/// Polygon is closed (last point -> next = first point)
/// 
/// https://doc.cgal.org/latest/Polyline_simplification_2/index.html
/// </summary>
/// <param name="threshold"></param>
/// <returns></returns>
Polygon Polygon::CreateSimplified(double threshold) const
{
	PolyLine simpleLine = PolyLine::CreateSimplified(threshold);
	if (simpleLine.Size() == 0)
	{
		return Polygon();
	}
	return Polygon(std::move(simpleLine));
}
