#ifndef MY_MATH_UTILS_H
#define MY_MATH_UTILS_H

#include <cstdint>
#include <type_traits>

#include "./Vector2.h"

#define MATH_RET_VAL(K) \
	typename std::enable_if<std::is_floating_point<K>::value, K>::type


#ifdef MYMATH_NAMESPACE
namespace MYMATH_NAMESPACE {
#endif
	struct MathUtils
	{		
		//==============================================================
		// Constants
		//==============================================================

		static const float EPSILON;
		static const float PI;
		static const uint32_t POWER_OF_TWO[32];
			
		//==============================================================
		// Classic methods
		//==============================================================
			
		static int FastPow10(int n);

		static int32_t FastRound(double val);
		static double RoundToDecimalPlaces(double val, int places);
		static double RoundTo(double val, double precision);
		static double RoundToNearest(double number, double multiple);

		static bool IsPowerOfTwo(uint64_t x);
		static bool IsZero(const float & value);
		static bool IsZero(const double & value);
				
		static float LineSegmentPointDistanceSquared(const Vector2f & a, const Vector2f & b, const Vector2f & p, float & t);
		static double LineSegmentPointDistanceSquared(const Vector2d & a, const Vector2d & b, const Vector2d & p, double & t);


		//==============================================================
		// Templated methods
		//==============================================================

		template <typename T>
		static bool IsInInterval(T start, T end, T val);

		template <typename T, typename K = T>
		static MATH_RET_VAL(K) DegToRad(const T & value);

		template <typename T, typename K = T>
		static MATH_RET_VAL(K) RadToDeg(const T & value);

		template <typename T, typename K = T>
		static K Signum(T val);

		template <typename T, typename K = T>
		static K Lerp(T start, T end, float value);

		template <typename T, typename K = T>
		static K Lerp(T start, T end, double value);

		template <typename T>
		static T BilinearInterpolation(
			const T &tx,
			const T &ty,
			const T &c00,
			const T &c10,
			const T &c01,
			const T &c11);
		
		template <typename T, typename K = T>
		static MATH_RET_VAL(K) CubicInterpolation(const T &A, const T &B, const T &C, const T &D, const T& pos);
		
		template <typename T, typename K = T>
		static MATH_RET_VAL(K) BiCubicInterpolation(const T(&values)[16], const T & x, const T & y);

		template <typename T, typename K = T>
		static MATH_RET_VAL(K) MapRange(const T &from1, const T &from2, const T &to1, const T &to2, const T &s);

		template <typename T, typename K = T>
		static MATH_RET_VAL(K) MapRangeFrom01(const T &to1, const T &to2, const T &s);

		template <typename T, typename K = T>
		static MATH_RET_VAL(K) MapRangeTo01(const T &from1, const T &from2, const T &s);
	};


	template <typename T>
	bool MathUtils::IsInInterval(T start, T end, T val)
	{
		return ((val >= start) && (val <= end));
	}

	template <typename T, typename K>
	MATH_RET_VAL(K) MathUtils::RadToDeg(const T & value)
	{
		return (value * K(57.2957795));
	};

	template <typename T, typename K>
	MATH_RET_VAL(K) MathUtils::DegToRad(const T & value)
	{
		return (value * K(0.0174532925));
	};

	template <typename T, typename K>
	K MathUtils::Signum(T value)
	{
		return static_cast<K>((T(0) < value) - (value < T(0)));
	};
	
	template <typename T>
	T MathUtils::BilinearInterpolation(
		const T &tx,
		const T &ty,
		const T &c00,
		const T &c10,
		const T &c01,
		const T &c11)
	{
#if 1
		//a = Lerp(c00, c10, tx)
		//b = Lerp(c01, c11, tx)
		//return Lerp(a, b, ty)
		T a = c00 * (1 - tx) + c10 * tx;
		T b = c01 * (1 - tx) + c11 * tx;
		return a * (1 - ty) + b * ty;
#else
		return (1 - tx) * (1 - ty) * c00 +
			tx * (1 - ty) * c10 +
			(1 - tx) * ty * c01 +
			tx * ty * c11;
#endif
	};
	
	template <typename T, typename K>
	K MathUtils::Lerp(T start, T end, float value)
	{
		return K((1.0f - value) * start + value * end);
	};

	template <typename T, typename K>
	K MathUtils::Lerp(T start, T end, double value)
	{
		return K((1.0 - value) * start + value * end);
	};


	/// <summary>
	/// http://www.paulinternet.nl/?page=bicubic
	/// </summary>
	/// <param name="A"></param>
	/// <param name="B"></param>
	/// <param name="C"></param>
	/// <param name="D"></param>
	/// <param name="pos"></param>
	/// <returns></returns>
	template <typename T, typename K>
	MATH_RET_VAL(K) MathUtils::CubicInterpolation(const T &A, const T &B, const T &C, const T &D, const T &pos)
	{
		return K(B + 0.5 * pos * (C - A + pos * (2.0 * A - 5.0 * B + 4.0 * C - D + pos * (3.0 * (B - C) + D - A))));
	};
	
	/// <summary>
	/// [x, y] must be in (0,1) interval
	/// </summary>
	/// <param name="values"></param>
	/// <param name="x"></param>
	/// <param name="y"></param>
	/// <returns></returns>
	template <typename T, typename K>
	MATH_RET_VAL(K) MathUtils::BiCubicInterpolation(const T(&values)[16], const T & x, const T & y)
	{
		K Ai = CubicInterpolation(values[0], values[1], values[2], values[3], x);
		K Bi = CubicInterpolation(values[4], values[5], values[6], values[7], x);
		K Ci = CubicInterpolation(values[8], values[9], values[10], values[11], x);
		K Di = CubicInterpolation(values[12], values[13], values[14], values[15], x);
		return CubicInterpolation(Ai, Bi, Ci, Di, y);
	};

	/// <summary>
	/// Map s from (from1, from2) range to (to1, to2)
	/// </summary>
	/// <param name="from1"></param>
	/// <param name="from2"></param>
	/// <param name="to1"></param>
	/// <param name="to2"></param>
	/// <param name="s"></param>
	/// <returns></returns>
	template <typename T, typename K>
	MATH_RET_VAL(K) MathUtils::MapRange(const T &from1, const T &from2, const T &to1, const T &to2, const T &s)
	{
		return to1 + (s - from1) * (to2 - to1) / K(from2 - from1);
	};

	/// <summary>
	/// Map s from (0, 1) range to (to1, to2)
	/// </summary>	
	/// <param name="to1"></param>
	/// <param name="to2"></param>
	/// <param name="s"></param>
	/// <returns></returns>
	template <typename T, typename K>
	MATH_RET_VAL(K) MathUtils::MapRangeFrom01(const T &to1, const T &to2, const T &s)
	{
		return to1 + (s) * (to2 - to1);
	};

	/// <summary>
	/// Map s from (from1, from2) range to (0, 1)
	/// </summary>	
	/// <param name="from1"></param>
	/// <param name="from2"></param>
	/// <param name="s"></param>
	/// <returns></returns>
	template <typename T, typename K>
	MATH_RET_VAL(K) MathUtils::MapRangeTo01(const T &from1, const T &from2, const T &s)
	{
		return (s - from1) * K(1) / K(from2 - from1);
	};

#undef MATH_RET_VAL

#ifdef MYMATH_NAMESPACE
}
#endif

#endif