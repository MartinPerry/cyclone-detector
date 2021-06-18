#ifndef SIMPLIFICATION_DOUGLAS_PEUCKER_H
#define SIMPLIFICATION_DOUGLAS_PEUCKER_H


#include <vector>

#include "./PolyLine.h"
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
		class SimplificationDouglasPeucker
		{
		public:
			SimplificationDouglasPeucker(const std::vector<MY_MATH_NS::Vector2d> & pts);

			PolyLine Simplify(double threshold) const;
		protected:

			const std::vector<MY_MATH_NS::Vector2d> & pts;

			double DistanceSquared(const MY_MATH_NS::Vector2d & a, const MY_MATH_NS::Vector2d & b) const;
			double PerpedincularDistanceSquared(const MY_MATH_NS::Vector2d & a, const MY_MATH_NS::Vector2d & b, const MY_MATH_NS::Vector2d & p) const;

			void Simplify(int startIndex, int endIndex, double tresshold, std::vector<int> & result) const;


		};
	}

#ifdef GRAPHICS_NAMESPACE
}
#endif

#undef MY_MATH_NS

#endif
