#include "PolyLine.h"

#include <cmath>

#include "./SimplificationDouglasPeucker.h"

#include "../../Math/MathUtils.h"

#ifdef MYMATH_NAMESPACE
using namespace MYMATH_NAMESPACE;
#endif


#ifdef GRAPHICS_NAMESPACE
using namespace GRAPHICS_NAMESPACE;
#endif

using namespace d2;

/// <summary>
/// ctor
/// </summary>
PolyLine::PolyLine()
{
}

/// <summary>
/// ctor from list of points
/// </summary>
/// <param name="pts"></param>
PolyLine::PolyLine(const std::vector<Vector2d> & pts)
{
	for (auto pt : pts)
	{
		this->AddPoint(pt);
	}
}

/// <summary>
/// copy ctor
/// </summary>
/// <param name="l"></param>
PolyLine::PolyLine(const PolyLine & l) : 
	bb(l.bb),
	pts(l.pts)
{	
}

/// <summary>
/// move ctor
/// </summary>
/// <param name="l"></param>
PolyLine::PolyLine(PolyLine && l) noexcept :
	bb(std::move(l.bb)),
	pts(std::move(l.pts))
{	
}

/// <summary>
/// Copy assignment
/// </summary>
/// <param name="l"></param>
/// <returns></returns>
PolyLine& PolyLine::operator=(const PolyLine& l)
{
	// Self-assignment detection
	if (&l == this)
	{
		return *this;
	}

	// Release any resource we're holding
	

	// Copy the resource
	this->bb = l.bb;
	this->pts = l.pts;

	return *this;
}

/// <summary>
/// Clear all points
/// </summary>
void PolyLine::Clear() noexcept
{
	this->bb = Aabb<Vector2d>();
	this->pts.clear();
}

/// <summary>
/// Insert point between existing points
/// </summary>
/// <param name="pt"></param>
/// <param name="atIndex"></param>
void PolyLine::InsertPoint(const Vector2d & pt, size_t atIndex)
{
	this->pts.insert(this->pts.begin() + atIndex, pt);
	bb.Update(pt);
}

/// <summary>
/// Erase single point at index
/// </summary>
/// <param name="atIndex"></param>
void PolyLine::ErasePoint(size_t atIndex)
{
	this->pts.erase(this->pts.begin() + atIndex);
}

/// <summary>
/// Add new point to the end line
/// </summary>
/// <param name="pt"></param>
void PolyLine::AddPoint(const Vector2d & pt)
{
	this->pts.push_back(pt);
	bb.Update(pt);	
}

/// <summary>
/// Add new point to the end of line if this point is different from the
/// last one added
/// </summary>
/// <param name="pt"></param>
/// <returns></returns>
bool PolyLine::AddPointNonDuplicate(const Vector2d & pt)
{
	if ((this->pts.size() > 0) && (this->pts.back() == pt))
	{
		//do not add duplicite points
		return false;
	}
	this->AddPoint(pt);
	return true;
}

/// <summary>
/// Remove last point from polyline
/// Bounding box is not recomputed !!!
/// </summary>
void PolyLine::RemoveLastPoint()
{
	this->pts.pop_back();
}

/// <summary>
/// Recompute bounding box
/// </summary>
void PolyLine::UpdateBoundingBox()
{
	this->bb = Aabb<Vector2d>();
	for (auto p : this->pts)
	{
		this->bb.Update(p);
	}
}

/// <summary>
/// Number of points in line
/// </summary>
/// <returns></returns>
size_t PolyLine::Size() const noexcept
{
	return this->pts.size();
}

/// <summary>
/// Computed Squared Lengh of polyline between point [from] and [to]
/// </summary>
/// <param name="from"></param>
/// <param name="to"></param>
/// <returns></returns>
double PolyLine::LengthSquaredFromTo(size_t from, size_t to) const
{
	if ((from == to) || (to == 0))
	{
		return 0;
	}
	
	if (from > to)
	{
		std::swap(to, from);
	}

	double len2 = 0;

	for (size_t i = from; i < to; i++)
	{
		double dx = this->pts[i + 1].x - this->pts[i].x;
		double dy = this->pts[i + 1].y - this->pts[i].y;


		len2 += ((dx * dx) + (dy * dy));
	}

	return len2;
}

/// <summary>
/// Computed Lengh of polyline between point [from] and [to]
/// </summary>
/// <param name="from"></param>
/// <param name="to"></param>
/// <returns></returns>
double PolyLine::LengthFromTo(size_t from, size_t to) const
{
	if ((from == to) || (to == 0))
	{
		return 0;
	}

	if (from > to)
	{
		std::swap(to, from);
	}

	double len2 = 0;

	for (size_t i = from; i < to; i++)
	{
		double dx = this->pts[i + 1].x - this->pts[i].x;
		double dy = this->pts[i + 1].y - this->pts[i].y;


		len2 += std::sqrt((dx * dx) + (dy * dy));
	}

	return len2;
}


/// <summary>
/// Computed Squared Lengh of polyline
/// </summary>
/// <returns></returns>
double PolyLine::LengthSquared() const
{
	if (pts.size() == 0)
	{
		return 0;
	}
	return this->LengthSquaredFromTo(0, this->pts.size() - 1);
}

/// <summary>
/// Computed length of polyline
/// </summary>
/// <returns></returns>
double PolyLine::Length() const
{
	if (pts.size() == 0)
	{
		return 0;
	}
	return this->LengthFromTo(0, this->pts.size() - 1);
}

/// <summary>
/// Calculate minimal distance of pt to polyline
/// segmentStart / segmentEnd - start and end point index of closest segment
/// </summary>
/// <param name="pt"></param>
/// <param name="segmentStart">(out)</param>
/// <param name="segmentEnd">(out)</param>
/// <returns></returns>
double PolyLine::GetMinDist(const Vector2d & pt, 
	size_t & segmentStart, size_t & segmentEnd, double & segmentT) const
{
	double d = std::numeric_limits<double>::max();
	segmentStart = size_t(-1);
	segmentEnd = size_t(-1);

	if (this->pts.size() < 2)
	{			
		segmentT = -1;
		return d;
	}

	double t;
	for (size_t i = 0; i < this->pts.size() - 1; i++)
	{				
		double curD = MathUtils::LineSegmentPointDistanceSquared(this->pts[i], this->pts[i + 1], pt, t);
		if (curD < d)
		{
			segmentT = t;
			d = curD;
			segmentStart = i;
			segmentEnd = i + 1;
		}
	}

	return std::sqrt(d);
}


const std::vector<Vector2d> & PolyLine::GetPoints() const noexcept
{
	return this->pts;
}

Vector2d PolyLine::GetPoint(size_t & segmentStart, double & segmentT) const
{
	auto dir = this->pts[segmentStart + 1] - this->pts[segmentStart];
	
	return this->pts[segmentStart] + dir * segmentT;
}

int PolyLine::PointOrientation(const Vector2d & pt) const
{
	size_t start, end;
	double t;
	this->GetMinDist(pt, start, end, t);

	double dx = this->pts[start].x - this->pts[end].x;
	double dy = this->pts[start].y - this->pts[end].y;

	double area = dx * (pt.y - this->pts[start].y) - dy * (pt.x - this->pts[start].x);

	int sign = 0;
	if (area < 0)
	{
		sign = -1;		
	}
	else if (area > 0)
	{
		sign = 1;
	}

	return sign;
/*
	double minArea = std::numeric_limits<double>::max();
	int minSign = 0;
	for (size_t i = 0; i < this->pts.size() - 1; i++)
	{
		double dx = this->pts[i].x - this->pts[i + 1].x;
		double dy = this->pts[i].y - this->pts[i + 1].y;

		double area = dx * (pt.y - this->pts[i].y) - dy * (pt.x - this->pts[i].x);
		int sign = 0;
		if (area < 0)
		{
			sign = -1;
			area *= -1;
		}
		else if (area > 0)
		{
			sign = 1;
		}

		if (area < minArea)
		{
			minArea = area;
			minSign = sign;
		}		
	}

	return minSign;
 */
}

/// <summary>
/// Scale line by (sx, sy) scale factors
/// </summary>
/// <param name="sx"></param>
/// <param name="sy"></param>
void PolyLine::Scale(double sx, double sy)
{
	this->bb.min.x *= sx;
	this->bb.min.y *= sy;

	this->bb.max.x *= sx;
	this->bb.max.y *= sy;

	for (size_t i = 0; i < this->pts.size(); i++)
	{
		this->pts[i].x *= sx;
		this->pts[i].y *= sy;
	}
}

/// <summary>
/// Test if any part of line is inside the box
/// </summary>
/// <param name="box"></param>
/// <returns></returns>
bool PolyLine::IsInBox(const Aabb<Vector2d> & box) const
{
	if (box.Intersect(this->bb) == false)
	{
		//line is entirely outside		
		return false;
	}

	if (box.IsInside(this->bb))
	{
		//line is entirely in	
		return true;
	}

	//line intersect bounding box

	size_t i = 0;
	Vector2d curPt = this->pts[i];
	int curCode = box.GetCode(curPt);

	do
	{
		Vector2d nextPt = this->pts[i + 1];		
		int nextCode = box.GetCode(nextPt);

		if (curCode & nextCode)
		{
			// bitwise AND is not 0: both points are outside window; we know nothing
			curPt = nextPt;
			curCode = nextCode;
			i++;
		}
		else
		{
			//intersection
			//or at least two points inside - accept
			return true;
		}		

	} while (i < this->pts.size() - 1);

	return false;
}

bool PolyLine::AddClampPoint(const Vector2d & pt)
{
	if ((this->pts.size() > 0) && (this->pts.back() == pt))
	{
		//do not add duplicite points
		return false;
	}
	this->pts.push_back(pt);
	return true;
}

/// <summary>
/// Clamp line using Cohen-Sutherland
/// </summary>
/// <param name="box"></param>
std::list<PolyLine> PolyLine::Clamp(const Aabb<Vector2d> & box) const
{
	std::vector<ClampInfo> ci;
	return this->ClampInternal(box, ci);
}

/// <summary>
/// Clamp line using Cohen-Sutherland
/// </summary>
/// <param name="box"></param>
/// <param name="clampInfo"></param>
/// <returns></returns>
std::list<PolyLine> PolyLine::ClampInternal(Aabb<Vector2d> box,
	std::vector<ClampInfo> & clampInfo) const
{
	clampInfo.clear();
	std::list<PolyLine> clampedLines;

	if (box.Intersect(this->bb) == false)
	{
		//line is entirely outside
		//return nothing
		return clampedLines;
	}
	
		
	//if the line is enclosed
	//end is set to be head->next
	//we need next, because first point is used as enclosed point
	//and must be processed
	//also, we need to use do - while
	//const SimpleLinePoint<Pixel> * const END_PTR = (head->prev != nullptr) ? head->next : nullptr;

	
	
	ClampInfo ci = { -1, -1 };
			
	if (this->bb.IsInside(box))
	{
		//line is entirely in		
		//no clipping needed

		clampedLines.push_back(*this);
		clampInfo.push_back(ci);	

		return clampedLines;
	}
	
	
	PolyLine cl;
	cl.bb = box;

	size_t it = 0;
	Vector2d curPt = this->pts[it];

	//line intersect bounding box
	//do the clipping
	do
	{
		Vector2d nextPt = this->pts[it + 1];

		ClampInfo curCi = { -1, -1 };
		if (PolyLine::ClampLineSegment(box, curPt, nextPt, curCi))
		{
			if ((curCi.in == -1) && (curCi.out == -1))
			{
				//both points were inside - no clamping
				//ad current
				cl.AddClampPoint(curPt);
			}
			else
			{
				//one/both points were outside - clamping

				if (curCi.in == -1)
				{
					cl.AddClampPoint(curPt);
					cl.AddClampPoint(nextPt);

					if (cl.pts.size() > 0)
					{
						//cur is inside, next was outside
						ci.out = curCi.out;

						clampInfo.push_back(ci);
						clampedLines.push_back(std::move(cl));
					}

					//empty line
					ci = { -1, -1 };
					cl = PolyLine();
					cl.bb = box;
				}
				else if (curCi.out == -1)
				{
					//cur was outside, next is inside
					ci.in = curCi.in;

					if (cl.pts.size() > 0)
					{
						clampInfo.push_back(ci);
						//last line is not empty					
						//add it to clamped lines
						clampedLines.push_back(std::move(cl));

						//empty line
						ci = { -1, -1 };
						cl = PolyLine();
						cl.bb = box;
					}


					//add clamped current point
					cl.AddClampPoint(curPt);
				}
				else
				{
					//both were outside - both are clamped
					//they form single line

					cl.AddClampPoint(curPt);
					cl.AddClampPoint(nextPt);

					if (cl.pts.size() > 0)
					{
						clampInfo.push_back(curCi);

						//end line
						clampedLines.push_back(std::move(cl));
					}

					//empty line
					ci = { -1, -1 };
					cl = PolyLine();
					cl.bb = box;
				}
			}
		}

		it++;
		curPt = this->pts[it];


	} while (it < this->Size() - 1);

	//add last line

	//add last point
	cl.AddClampPoint(curPt);
	if (cl.pts.size() > 0)
	{
		//line has at least two points = proper line					
		clampInfo.push_back(ci);
		clampedLines.push_back(std::move(cl));
	}
		
	//----------------------------------------------------------------------
	//check first and last segment
	//this is more-or-less used only if line is forming enclosed polygon
	//check only if there is more than one line
	//because if there is one line, first and last point may be the same -> if line is totaly in and enclosed
	if ((clampedLines.size() > 1) && (clampedLines.front().pts.front() == clampedLines.back().pts.back()))
	{
		//first and last segments are not connected

		auto & firstSegment = clampedLines.front();
		size_t firstSegmentPtsCount = firstSegment.pts.size();

		//copy points from first segment to last segment
		//first point is not copied, since 
		//last segment last point = first segment first point
		for (size_t i = 1; i < firstSegmentPtsCount; i++)
		{
			clampedLines.back().pts.push_back(firstSegment[i]);
		}

		//set clamping info for last segment
		//output is same as output of first segment
		clampInfo.back().out = clampInfo[0].out;

		//remove first segment
		clampedLines.pop_front();	
		clampInfo.erase(clampInfo.begin());
	}
	//--------------------------------------------------------------------	
	
	//remove degenarted "lines" that contains only one point
	size_t clampIndex = 0;
	std::list<PolyLine>::iterator cit = clampedLines.begin();
	while (cit != clampedLines.end())
	{		
		if (cit->Size() == 1)
		{
			clampInfo.erase(clampInfo.begin() + clampIndex);
			cit = clampedLines.erase(cit);
		}
		else
		{
			clampIndex++;
			cit++;
		}		
	}
	
	return clampedLines;
}

/// <summary>
/// Clamp line segment curPt - nextPt
/// and update its points
/// ClampInfo contains info about entry / exit edge after clamping
/// If points after clamping are not in bounding box, method return false
/// Using Cohen-Sutherland algorithm from:
/// https://en.wikipedia.org/wiki/Cohen%E2%80%93Sutherland_algorithm
/// </summary>
/// <param name="box"></param>
/// <param name="curPt"></param>
/// <param name="nextPt"></param>
/// <param name="ci"></param>
/// <returns></returns>
bool PolyLine::ClampLineSegment(const Aabb<Vector2d> & box,
	Vector2d & curPt, Vector2d & nextPt, ClampInfo & ci)
{
	int edgeInfo = -1;
	int curCode = box.GetCode(curPt);
	int nextCode = box.GetCode(nextPt);

	while (true)
	{
		if (!(curCode | nextCode))
		{
			// bitwise OR is 0: both points inside window; trivially accept and exit loop
			return true;
		}
		else if (curCode & nextCode)
		{
			// bitwise AND is not 0: both points share an outside zone (LEFT, RIGHT, TOP,
			// or BOTTOM), so both must be outside window; exit loop (accept is false)
			return false;
		}
		else
		{
			// failed both tests, so calculate the line segment to clip
			// from an outside point to an intersection with clip edge
			Vector2d p;

			// At least one endpoint is outside the clip rectangle; pick it.
			int codeOut = curCode ? curCode : nextCode;

			if (codeOut & Aabb<Vector2d>::TOP)
			{
				// point is above the clip window
				p.x = curPt.x + (nextPt.x - curPt.x) * (box.max.y - curPt.y) / (nextPt.y - curPt.y);
				p.y = box.max.y;
				edgeInfo = Aabb<Vector2d>::CLAMP_TOP;
			}
			else if (codeOut & Aabb<Vector2d>::BOTTOM)
			{
				// point is below the clip window
				p.x = curPt.x + (nextPt.x - curPt.x) * (box.min.y - curPt.y) / (nextPt.y - curPt.y);
				p.y = box.min.y;
				edgeInfo = Aabb<Vector2d>::CLAMP_BOTTOM;
			}
			else if (codeOut & Aabb<Vector2d>::RIGHT)
			{
				// point is to the right of clip window
				p.y = curPt.y + (nextPt.y - curPt.y) * (box.max.x - curPt.x) / (nextPt.x - curPt.x);
				p.x = box.max.x;
				edgeInfo = Aabb<Vector2d>::CLAMP_RIGHT;
			}
			else if (codeOut & Aabb<Vector2d>::LEFT)
			{
				// point is to the left of clip window
				p.y = curPt.y + (nextPt.y - curPt.y) * (box.min.x - curPt.x) / (nextPt.x - curPt.x);
				p.x = box.min.x;
				edgeInfo = Aabb<Vector2d>::CLAMP_LEFT;
			}

			if (codeOut == curCode)
			{
				curPt = p;
				ci.in = edgeInfo;
				curCode = box.GetCode(curPt);
			}
			else
			{
				nextPt = p;
				ci.out = edgeInfo;
				nextCode = box.GetCode(nextPt);
			}
		}
	}
}


/// <summary>
/// Simplify line
/// </summary>
/// <param name="threshold"></param>
/// <returns></returns>
PolyLine PolyLine::CreateSimplified(double threshold) const
{	
	SimplificationDouglasPeucker s(this->pts);
	return s.Simplify(threshold);
}
