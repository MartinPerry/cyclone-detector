#include "./MathUtils.h"

#include <cmath>

#ifdef MYMATH_NAMESPACE
using namespace MYMATH_NAMESPACE;
#endif

const float MathUtils::EPSILON = 0.00001f;
const float MathUtils::PI = 3.1415927410125732421875f;
const uint32_t MathUtils::POWER_OF_TWO[32] = {
	static_cast<uint32_t>(std::pow(2, 0)), static_cast<uint32_t>(std::pow(2, 1)),
	static_cast<uint32_t>(std::pow(2, 2)), static_cast<uint32_t>(std::pow(2, 3)),
	static_cast<uint32_t>(std::pow(2, 4)), static_cast<uint32_t>(std::pow(2, 5)),
	static_cast<uint32_t>(std::pow(2, 6)), static_cast<uint32_t>(std::pow(2, 7)),
	static_cast<uint32_t>(std::pow(2, 8)), static_cast<uint32_t>(std::pow(2, 9)),
	static_cast<uint32_t>(std::pow(2, 10)), static_cast<uint32_t>(std::pow(2, 11)),
	static_cast<uint32_t>(std::pow(2, 12)), static_cast<uint32_t>(std::pow(2, 13)),
	static_cast<uint32_t>(std::pow(2, 14)), static_cast<uint32_t>(std::pow(2, 15)),
	static_cast<uint32_t>(std::pow(2, 16)), static_cast<uint32_t>(std::pow(2, 17)),
	static_cast<uint32_t>(std::pow(2, 18)), static_cast<uint32_t>(std::pow(2, 19)),
	static_cast<uint32_t>(std::pow(2, 20)), static_cast<uint32_t>(std::pow(2, 21)),
	static_cast<uint32_t>(std::pow(2, 22)), static_cast<uint32_t>(std::pow(2, 23)),
	static_cast<uint32_t>(std::pow(2, 24)), static_cast<uint32_t>(std::pow(2, 25)),
	static_cast<uint32_t>(std::pow(2, 26)), static_cast<uint32_t>(std::pow(2, 27)),
	static_cast<uint32_t>(std::pow(2, 28)), static_cast<uint32_t>(std::pow(2, 29)),
	static_cast<uint32_t>(std::pow(2, 30)), static_cast<uint32_t>(std::pow(2, 31))
};


//==========================================================================
// Rounding
//==========================================================================

/// <summary>
/// Calculate 10^n
/// https://stackoverflow.com/questions/18581560/any-way-faster-than-pow-to-compute-an-integer-power-of-10-in-c
/// </summary>
/// <param name="n"></param>
/// <returns></returns>
int MathUtils::FastPow10(int n)
{
	int r = 1;
	const int x = 10;
	while (n)
	{
		if (n & 1)
		{
			r *= x;
			n--;
		}
		else
		{
			r *= x * x;
			n -= 2;
		}
	}

	return r;
}

/// <summary>
/// Fast round of double to 32bit int
/// 0.5 is rounded to nearest int (4.5 => 4; 4.51 => 5)
/// BEWARE - round it to 32bit int !!!! precision is limited
/// Source: https://stackoverflow.com/questions/17035464/a-fast-method-to-round-a-double-to-a-32-bit-int-explained
/// </summary>
/// <param name="val"></param>
/// <returns></returns>
int32_t MathUtils::FastRound(double val)
{
	val += 6755399441055744.0;
	return reinterpret_cast<int32_t&>(val);
}

/// <summary>
/// Round number to a given precision (number of decimal places)
/// </summary>
/// <param name="val"></param>
/// <param name="precision"></param>
/// <returns></returns>
double MathUtils::RoundToDecimalPlaces(double val, int places)
{	
	double precision = 1.0 / FastPow10(places);
	return MathUtils::RoundTo(val, precision);
}

/// <summary>
/// Round number to a given precision
/// Precision is a decimal number - 
/// > 0 - rounding in whole part
/// < 0 - rounding in fractional part
/// </summary>
/// <param name="val"></param>
/// <param name="precision"></param>
/// <returns></returns>
double MathUtils::RoundTo(double val, double precision)
{	
	val /= precision;
	return std::floor(val + 0.5) * precision;
}

/// <summary>
/// Round number to nearest multiple
/// Eg: 
/// input: 0.3, multiple: 0.2 => returns 0.4
/// input: 0.45, multiple: 0.2 => returns 0.4
/// </summary>
/// <param name="number"></param>
/// <param name="multiple"></param>
/// <returns></returns>
double MathUtils::RoundToNearest(double number, double multiple)
{
	double half = multiple / 2;
	return number + half - std::fmod(number + half, multiple);
}


//==========================================================================
// Tests
//==========================================================================

/// <summary>
/// Test if number is power of two
/// </summary>
/// <param name="x"></param>
/// <returns></returns>
bool MathUtils::IsPowerOfTwo(uint64_t x)
{
	/* While x is even and > 1 */
	while (((x % 2) == 0) && x > 1)
	{
		x /= 2;
	}
	return (x == 1);
};

bool MathUtils::IsZero(const float & value)
{
	if (value < -EPSILON) return false;
	if (value > EPSILON) return false;
	return true;
}

bool MathUtils::IsZero(const double & value)
{
	if (value < -EPSILON) return false;
	if (value > EPSILON) return false;
	return true;
}


float MathUtils::LineSegmentPointDistanceSquared(const Vector2f & a, const Vector2f & b, const Vector2f & p, float & t)
{
	// Return minimum distance between line segment vw and point p
	const float l2 = Vector2f::DistanceSquared(a, b);  // i.e. |w-v|^2 -  avoid a sqrt
	if (IsZero(l2))
	{
		t = 0.0f;
		return Vector2f::DistanceSquared(p, a);   // v == w case
	}

	// Consider the line extending the segment, parameterized as v + t (w - v).
	// We find projection of point p onto the line. 
	// It falls where t = [(p-v) . (w-v)] / |w-v|^2
	t = Vector2f::Dot(p - a, b - a) / l2;
	if (t < 0.0f)
	{
		return Vector2f::DistanceSquared(p, a);       // Beyond the 'v' end of the segment
	}
	else if (t > 1.0f)
	{
		return Vector2f::DistanceSquared(p, b);  // Beyond the 'w' end of the segment
	}

	Vector2f projection = a + t * (b - a);  // Projection falls on the segment
	return Vector2f::DistanceSquared(p, projection);
}

double MathUtils::LineSegmentPointDistanceSquared(const Vector2d & a, const Vector2d & b, const Vector2d & p, double & t)
{
	// Return minimum distance between line segment vw and point p
	const double l2 = Vector2d::DistanceSquared(a, b);  // i.e. |w-v|^2 -  avoid a sqrt
	if (IsZero(l2))
	{
		t = 0.0;
		return Vector2d::DistanceSquared(p, a);   // v == w case
	}

	// Consider the line extending the segment, parameterized as v + t (w - v).
	// We find projection of point p onto the line. 
	// It falls where t = [(p-v) . (w-v)] / |w-v|^2
	t = Vector2d::Dot(p - a, b - a) / l2;
	if (t < 0.0)
	{
		return Vector2d::DistanceSquared(p, a);       // Beyond the 'v' end of the segment
	}
	else if (t > 1.0)
	{
		return Vector2d::DistanceSquared(p, b);  // Beyond the 'w' end of the segment
	}

	Vector2d projection = a + t * (b - a);  // Projection falls on the segment
	return Vector2d::DistanceSquared(p, projection);
}