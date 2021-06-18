#ifndef VECTOR2_H
#define VECTOR2_H


#include <limits>

#ifdef MYMATH_NAMESPACE
namespace MYMATH_NAMESPACE {
#endif
	
	template <typename T = float>
	struct Vector2
	{
	public:

		union {
			struct {
				T x;
				T y;
			};
			T v[2];
		};

		//ctors
		Vector2() noexcept;
		Vector2(float x, float y) noexcept;
		Vector2(double x, double y) noexcept;
		Vector2(int x, int y) noexcept;
		Vector2(const Vector2<T> & other) noexcept;

		//getters
		//T* ToArray() {  T p[3] = {X, Y, Z}; return p; };

		static Vector2<T> UnitX() noexcept { return Vector2<T>(1, 0); };
		static Vector2<T> UnitY() noexcept { return Vector2<T>(0, 1); };
		static Vector2<T> Zero() noexcept { return Vector2<T>(0, 0); };
		static Vector2<T> One() noexcept { return Vector2<T>(1, 1); };
		static Vector2<T> Infinity() noexcept { return Vector2<T>(std::numeric_limits<T>::max(), std::numeric_limits<T>::max()); };
		static int SizeInBytes() noexcept { return 2 * sizeof(T); };
		static int ElementsCount() noexcept { return 2; };

		//static functions			
		static T CrossSquared(const Vector2<T> &left, const Vector2<T> &right) noexcept;
		static Vector2<T> Center(const Vector2<T> &left, const Vector2<T> &right) noexcept;
		static T Dot(const Vector2<T> &left, const Vector2<T> &right) noexcept;
		static T Distance(const Vector2<T> &first, const Vector2<T> &second) noexcept;
		static T DistanceSquared(const Vector2<T> &first, const Vector2<T> &second) noexcept;		
		static Vector2<T> Normalize(const Vector2<T> &vector) noexcept;
		static Vector2<T> MakeOrtogonal(const Vector2<T> & v) noexcept;

		static int SortCompare(const Vector2<T> & p1, const Vector2<T> & p2) noexcept;
		
		//functions
		T Length() const noexcept;
		T LengthSquared() const noexcept;
		T Distance(const Vector2<T> &value) const noexcept;
		T DistanceSquared(const Vector2<T> &value) const noexcept;
		void Normalize() noexcept;

		bool IsZero() const noexcept;
		bool IsInfinity() const noexcept;

		Vector2<T> operator +(const Vector2<T> &value) const noexcept;
		Vector2<T> operator -(const Vector2<T> &value) const noexcept;
		Vector2<T> operator -() const noexcept;
		Vector2<T> operator *(T value) const noexcept;
		Vector2<T> operator /(T value) const noexcept;
		T & operator [](const int index);

		void operator +=(const Vector2<T> &value) noexcept;
		void operator -=(const Vector2<T> &value) noexcept;
		void operator *=(T value) noexcept;

		bool operator ==(const Vector2<T> &value) const noexcept;
		bool operator !=(const Vector2<T> &value) const noexcept;



		friend Vector2<T> operator /(T scale, const Vector2<T> &value) noexcept
		{
			return Vector2<T>(scale / value.x, scale / value.y);
		}

		friend Vector2<T> operator *(T value, const Vector2<T> &right) noexcept
		{
			return Vector2<T>(value * right.x, value * right.y);
		}			
	};


	template <typename T>
	inline Vector2<T> Vector2<T>::operator +(const Vector2<T> &value) const noexcept
	{
		return Vector2<T>(x + value.x, y + value.y);
	}

	template <typename T>
	inline Vector2<T> Vector2<T>::operator -(const Vector2<T> &value) const noexcept
	{
		return Vector2(x - value.x, y - value.y);
	}

	template <typename T>
	inline Vector2<T> Vector2<T>::operator -() const noexcept
	{
		return Vector2<T>(-x, -y);
	}

	template <typename T>
	inline Vector2<T> Vector2<T>::operator /(T value) const noexcept
	{
		T r = T(1.0) / value;
		return Vector2<T>(x * r, y * r);
	}

	template <typename T>
	inline Vector2<T> Vector2<T>::operator *(T value) const noexcept
	{
		return Vector2<T>(x * value, y * value);
	}

	template <typename T>
	inline void Vector2<T>::operator +=(const Vector2 &value) noexcept
	{
		x += value.x;
		y += value.y;
	}

	template <typename T>
	inline void Vector2<T>::operator -=(const Vector2 &value) noexcept
	{
		x -= value.x;
		y -= value.y;
	}

	template <typename T>
	inline void Vector2<T>::operator *=(T value) noexcept
	{
		x *= value;
		y *= value;
	}


	template <typename T>
	inline bool Vector2<T>::operator ==(const Vector2 &value) const noexcept
	{
		return (x == value.x && y == value.y);
	}

	template <typename T>
	inline bool Vector2<T>::operator !=(const Vector2 &value) const noexcept
	{
		return (x != value.x || y != value.y);
	}

	template <typename T>
	inline T & Vector2<T>::operator [](const int index)
	{
		return v[index];
	}


	using Vector2f = Vector2<float>;
	using Vector2d = Vector2<double>;

#ifdef MYMATH_NAMESPACE
};
#endif

//=========================================================================================================
//=========================================================================================================
//=========================================================================================================

#include <unordered_map>

#ifdef MYMATH_NAMESPACE
using namespace MYMATH_NAMESPACE;
#endif

/// <summary>
/// For use in std::unordered_map
/// http://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
/// </summary>
namespace std
{
	template <>
	struct hash<Vector2<float>>
	{
		std::size_t operator()(const Vector2<float> & k) const
		{
			//Simple hash - convert float to int bit representation
			union Tmp
			{
				float valFloat;
				int32_t valInt32;
			};

			Tmp x, y;

			x.valFloat = k.x;
			y.valFloat = k.y;
			
			std::size_t hash = 17;
			hash = hash * 23 + x.valInt32;
			hash = hash * 23 + y.valInt32;
			return hash;
		};
	};

	template <>
	struct hash<Vector2<double>>
	{
		std::size_t operator()(const Vector2<double> & k) const
		{
			//Simple hash - convert float to int bit representation
			union Tmp
			{
				double valDouble;
				int64_t valInt64;
			};

			Tmp x, y;

			x.valDouble = k.x;
			y.valDouble = k.y;

			int64_t hash = 17;
			hash = hash * 23 + x.valInt64;
			hash = hash * 23 + y.valInt64;
			return static_cast<std::size_t>(hash);
		};
	};
};

#endif