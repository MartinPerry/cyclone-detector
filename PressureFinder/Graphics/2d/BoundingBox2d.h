#ifndef MY_BOUNDING_BOX_2D_H
#define MY_BOUNDING_BOX_2D_H

#include <limits>
#include <optional>

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

		/// <summary>
		/// Axis-aligned bounding box in 2D
		/// </summary>
		template <typename T>
		struct Aabb
		{				
			static const int INSIDE = 0; // 0000
			static const int LEFT = 1;   // 0001
			static const int RIGHT = 2;  // 0010
			static const int BOTTOM = 4; // 0100
			static const int TOP = 8;    // 1000

			static const int CLAMP_TOP = 0;
			static const int CLAMP_LEFT = 1;
			static const int CLAMP_BOTTOM = 2;
			static const int CLAMP_RIGHT = 3;


			T min;
			T max;

			using Type = decltype(min.x);

			/// <summary>
			/// Create union of two AABBs
			/// </summary>
			/// <param name="a"></param>
			/// <param name="b"></param>
			/// <returns></returns>
			static Aabb CreateUnion(const Aabb<T> & a, const Aabb<T> & b)
			{				
				Type minX = std::min(a.min.x, b.min.x);
				Type minY = std::min(a.min.y, b.min.y);
				
				Type maxX = std::max(a.max.x, b.max.x);
				Type maxY = std::max(a.max.y, b.max.y);

				return Aabb<T>({ minX, minY }, { maxX, maxY });
			}

			/// <summary>
			/// Create union of AABB and single point
			/// </summary>
			/// <param name="a"></param>
			/// <param name="b"></param>
			/// <returns></returns>
			static Aabb CreateUnion(const Aabb<T> & a, const T & b)
			{
				Type minX = std::min(a.min.x, b.x);
				Type minY = std::min(a.min.y, b.y);

				Type maxX = std::max(a.max.x, b.x);
				Type maxY = std::max(a.max.y, b.y);

				return Aabb<T>({ minX, minY }, { maxX, maxY });
			}

			/// <summary>
			/// ctor
			/// Init min and max to default values
			/// resulting in an "non existing" box
			/// min => max
			/// max => min
			/// </summary>
			Aabb() noexcept
			{
				this->Clear();
			};

			/// <summary>
			/// ctor
			/// </summary>
			/// <param name="min"></param>
			/// <param name="max"></param>
			Aabb(const T & min, const T & max) noexcept :
				min(min),
				max(max)
			{
			};

			/// <summary>
			/// Clear and set min/max to default values
			/// resulting in an "non existing" box
			/// min => max
			/// max => min
			/// </summary>
			void Clear() noexcept
			{
				min.x = std::numeric_limits<Type>::max();
				min.y = std::numeric_limits<Type>::max();
				max.x = std::numeric_limits<Type>::lowest();
				max.y = std::numeric_limits<Type>::lowest();
			};

			/// <summary>
			/// Update bounding box with input point
			/// </summary>
			/// <param name="pt"></param>
			void Update(const T & pt) noexcept
			{
				if (min.x > pt.x) min.x = pt.x;
				if (min.y > pt.y) min.y = pt.y;
				if (max.x < pt.x) max.x = pt.x;
				if (max.y < pt.y) max.y = pt.y;
			};

			/// <summary>
			/// Enlarge box by percents (0 - 1 range) of each extens size
			/// </summary>
			/// <param name="percents"></param>
			void Enlarge(float percents) noexcept
			{
				Type dx = percents * this->GetWidth();
				Type dy = percents * this->GetHeight();

				min.x -= dx;
				max.x += dx;

				min.y -= dy;
				max.y += dy;
			};

			/// <summary>
			/// Get area of AABB
			/// </summary>
			/// <returns></returns>
			Type GetArea() const noexcept
			{
				return this->GetWidth() * this->GetHeight();
			};

			
			Type GetWidth() const noexcept
			{
				return max.x - min.x;
			};

			Type GetHeight() const noexcept
			{
				return max.y - min.y;
			};

			/// <summary>
			/// Get centroid of AABB
			/// </summary>
			/// <returns></returns>
			T GetCentroid() const noexcept
			{
				return { 0.5 * (max.x + min.x), 0.5 * (max.y + min.y) };
			}

			/// <summary>
			/// Get corner based on its ID of edge
			/// ID corresponds to CLAMP_xxx constant
			/// which is starting point of edge in CW order
			/// </summary>
			/// <param name="index"></param>
			/// <returns></returns>
			T GetCorner(int index) const noexcept
			{
				switch (index)
				{
				case CLAMP_TOP: return { max.x, max.y };
				case CLAMP_LEFT: return { min.x, max.y };
				case CLAMP_BOTTOM: return { min.x, min.y };
				default: return { max.x, min.y };
				}
			};

			/// <summary>
			/// Get maximal axis
			/// 0 - width is max
			/// 1 - height is max
			/// </summary>
			/// <returns></returns>
			int GetMaxExtent() const noexcept
			{				
				if (this->GetWidth() > this->GetHeight())
				{
					return 0;
				}
				else
				{
					return 1;
				}
			}

			/// <summary>
			/// Get area of intersection with bb as a new AABB
			/// If there is no intersection = return nullopt
			/// If one is fully inside other, return area of smaller
			/// </summary>
			/// <param name="bb"></param>
			/// <returns></returns>
			std::optional<Aabb<T>> GetIntersect(const Aabb<T> & bb) const noexcept
			{
				if (this->Intersect(bb) == false)
				{
					return std::nullopt;
				}
								
				Type minX = std::max(this->min.x, bb.min.x);
				Type minY = std::max(this->min.y, bb.min.y);
				
				Type maxX = std::min(this->max.x, bb.max.x);
				Type maxY = std::min(this->max.y, bb.max.y);

				return Aabb<T>({ minX, minY }, { maxX, maxY });
			}
			

			/// <summary>
			/// Check if bounding box "bb" intersects "this"
			/// </summary>
			/// <param name="bb"></param>
			/// <returns></returns>
			bool Intersect(const Aabb<T> & bb) const noexcept
			{
				if (bb.min.x > max.x) return false;
				if (bb.min.y > max.y) return false;

				if (bb.max.x < min.x) return false;
				if (bb.max.y < min.y) return false;

				return true;
			};

			/// <summary>
			/// Intersection with ray given by two points A, B
			/// https://tavianator.com/fast-branchless-raybounding-box-intersections/
			/// </summary>
			/// <param name="a"></param>
			/// <param name="b"></param>
			/// <returns></returns>
			bool IntersectRay(const T & a, const T & b) const  noexcept
			{
				double rayInvDirX = 1.0 / (b.x - a.x);
				double rayInvDirY = 1.0 / (b.y - a.y);

				double tx1 = (min.x - a.x) * rayInvDirX;
				double tx2 = (max.x - a.x) * rayInvDirX;

				double tmin = std::min(tx1, tx2);
				double tmax = std::max(tx1, tx2);

				double ty1 = (min.y - a.y) * rayInvDirY;
				double ty2 = (max.y - a.y) * rayInvDirY;

				tmin = std::max(tmin, std::min(ty1, ty2));
				tmax = std::min(tmax, std::max(ty1, ty2));

				return tmax >= tmin;
			};
			
			/// <summary>
			/// Chek if "this" bounding box is inside "bb"
			/// </summary>
			/// <param name="bb"></param>
			/// <returns></returns>
			bool IsInside(const Aabb<T> & bb) const noexcept
			{
				return ((min.x > bb.min.x) && (max.x < bb.max.x) &&
					(min.y > bb.min.y) && (max.y < bb.max.y));
			};

			/// <summary>
			/// Chek if point "p" is inside "this"
			/// </summary>
			/// <param name="p"></param>
			/// <returns></returns>
			bool IsInside(const T & p) const noexcept
			{
				return ((p.x > min.x) && (p.x < max.x) &&
					(p.y > min.y) && (p.y < max.y));
			};

			/// <summary>
			/// Get clamping code for Cohen-Sutherland
			/// </summary>
			/// <param name="p"></param>
			/// <returns></returns>
			int GetCode(const T & p) const noexcept
			{
				int code = INSIDE;

				if (p.x < min.x)   code |= LEFT;
				else if (p.x > max.x)  code |= RIGHT;

				if (p.y < min.y) code |= BOTTOM;
				else if (p.y > max.y) code |= TOP;

				return code;
			};
		};
	}
	
#ifdef GRAPHICS_NAMESPACE
}
#endif

#undef MY_MATH_NS

#endif