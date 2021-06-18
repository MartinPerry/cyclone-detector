#ifndef MY_POLYLINE_H
#define MY_POLYLINE_H


#include <vector>
#include <list>

#include "./BoundingBox2d.h"
#include "../../Math/Vector2.h"

#ifdef MYMATH_NAMESPACE 
#	define MY_MATH_NS MYMATH_NAMESPACE 
#else
#	define MY_MATH_NS 
#endif

#ifdef GRAPHICS_NAMESPACE
namespace GRAPHICS_NAMESPACE {
#endif

	namespace d2
	{
		class PolyLine
		{
		public:
			PolyLine();
			PolyLine(const std::vector<MY_MATH_NS::Vector2d> & pts);
			PolyLine(const PolyLine & l);
			PolyLine(PolyLine && l) noexcept;

			PolyLine& operator=(const PolyLine& a);

			void Clear() noexcept;

			void InsertPoint(const MY_MATH_NS::Vector2d & pt, size_t atIndex);
			void ErasePoint(size_t atIndex);

			void AddPoint(const MY_MATH_NS::Vector2d & pt);
			bool AddPointNonDuplicate(const MY_MATH_NS::Vector2d & pt);
			size_t Size() const noexcept;
			double Length() const;
			double LengthSquared() const;
			double LengthFromTo(size_t from, size_t to) const;
			double LengthSquaredFromTo(size_t from, size_t to) const;
			double GetMinDist(const MY_MATH_NS::Vector2d & pt, 
				size_t & segmentStart, size_t & segmentEnd, double & segmentT) const;

			const std::vector<MY_MATH_NS::Vector2d> & GetPoints() const noexcept;
			MY_MATH_NS::Vector2d GetPoint(size_t & segmentStart, double & segmentT) const;

			int PointOrientation(const MY_MATH_NS::Vector2d & pt) const;

			void Scale(double sx, double sy);

			void RemoveLastPoint();
			void UpdateBoundingBox();

			const MY_MATH_NS::Vector2d & GetLastPoint() const
			{
				return this->pts.back();
			}

			const Aabb<MY_MATH_NS::Vector2d> & GetBoundingBox() const noexcept
			{
				return this->bb;
			}

			MY_MATH_NS::Vector2d & operator[](size_t index)
			{
				return pts[index];
			};

			const MY_MATH_NS::Vector2d & operator[](size_t index) const
			{
				return pts[index];
			};

			bool IsInBox(const Aabb<MY_MATH_NS::Vector2d> & box) const;

			std::list<PolyLine> Clamp(const Aabb<MY_MATH_NS::Vector2d> & box) const;

			PolyLine CreateSimplified(double threshold) const;

			friend class Polygon;

		protected:

			struct ClampInfo
			{
				int in;
				int out;
			};

			//minimal bounding box. It is not quaranteed to be tight
			//eg. after clipping, bb is set to clip window
			Aabb<MY_MATH_NS::Vector2d> bb;

			mutable std::vector<MY_MATH_NS::Vector2d> pts;

			bool AddClampPoint(const MY_MATH_NS::Vector2d & pt);

			std::list<PolyLine> ClampInternal(Aabb<MY_MATH_NS::Vector2d> box,
				std::vector<ClampInfo> & clampInfo) const;

			static bool ClampLineSegment(const Aabb<MY_MATH_NS::Vector2d> & box,
				MY_MATH_NS::Vector2d & curPt, MY_MATH_NS::Vector2d & nextPt, ClampInfo & ci);


		};
	}

#ifdef GRAPHICS_NAMESPACE
}
#endif

#undef MY_MATH_NS

#endif