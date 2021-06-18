#ifndef MY_POLYGON_H
#define MY_POLYGON_H

#ifdef MYMATH_NAMESPACE
namespace MYMATH_NAMESPACE {
#endif	
	template <typename T>
	struct Vector2;
#ifdef MYMATH_NAMESPACE
}
#endif

#include <list>

#include "./PolyLine.h"

#ifdef MYMATH_NAMESPACE 
#	define MY_MATH_NS MYMATH_NAMESPACE 
#else
#	define MY_MATH_NS 
#endif

using Vector2d = MY_MATH_NS::Vector2<double>;

#ifdef GRAPHICS_NAMESPACE
namespace GRAPHICS_NAMESPACE {
#endif

	namespace d2
	{

		/// <summary>
		/// Polygon
		/// 
		/// Inherits from Line:
		/// Protected inheritance - we hide parent members Clamp and Simplify
		/// Without this, we could do Line l = Polygon()
		/// and l.Clamp() would call Lines Clamp method
		/// With protected inheritance, we disable this and 
		/// can no longer do Line l = Polygon()
		/// Also we dont need virtual methods, we just hide parent implementation
		/// which is faster
		/// </summary>
		class Polygon : virtual protected PolyLine
		{
		public:
			using PolyLine::Size;
			using PolyLine::InsertPoint;
			using PolyLine::AddPoint;
			using PolyLine::AddPointNonDuplicate;
			using PolyLine::GetBoundingBox;
			using PolyLine::operator[];
			using PolyLine::LengthFromTo;
			using PolyLine::LengthSquaredFromTo;

			static Polygon CreateFromPolyline(const PolyLine & l);

			Polygon() noexcept;
			Polygon(const Polygon & p);
			Polygon(Polygon && p) noexcept;
			Polygon& operator=(const Polygon& p);

			double CalcArea() const;

			bool IsPointInside(const Vector2d & p) const;

			Vector2d CalcCentroid() const;

			bool IsInBox(const Aabb<Vector2d> & box) const;
			std::list<Polygon> Clamp(const Aabb<Vector2d> & box)  const;


			Polygon CreateSimplified(double threshold) const;

		protected:

			struct ClampPixel : public Vector2d
			{
				ClampPixel(const Vector2d & d, bool corner) : Vector2d(d),
					corner(corner)
				{}

				bool corner;
			};

			Polygon(PolyLine && l);

			std::list<ClampPixel> BuildClipWindow(const std::list<PolyLine> & clampedLines,
				const std::vector<ClampInfo> & clampInfo, const Aabb<Vector2d> & box) const;

		};
	}

#ifdef GRAPHICS_NAMESPACE
}
#endif

#undef MY_MATH_NS

#endif