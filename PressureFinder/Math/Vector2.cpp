#include "./Vector2.h"

#include <cmath>


#ifdef MYMATH_NAMESPACE
using namespace MYMATH_NAMESPACE;
#endif

template <typename T>
Vector2<T>::Vector2(float x, float y) noexcept :
	x(x),
	y(y)
{		
}

template <typename T>
Vector2<T>::Vector2(double x, double y) noexcept :
	x(T(x)),
	y(T(y))
{
}
	
template <typename T>
Vector2<T>::Vector2(int x, int y) noexcept :
	x(T(x)),
	y(T(y))
{
}

template <typename T>
Vector2<T>::Vector2() noexcept :
	Vector2<T>(T(0.0), T(0.0))
{		
}

template <typename T>
Vector2<T>::Vector2(const Vector2<T> & other) noexcept :
	Vector2<T>(other.x, other.y)
{		
}


template <typename T>
T Vector2<T>::Length() const noexcept
{
	return std::sqrt( LengthSquared() );
}

template <typename T>
T Vector2<T>::LengthSquared() const noexcept
{
	return x * x + y * y;
}

template <typename T>
T Vector2<T>::Distance(const Vector2<T> &value) const noexcept
{	
	return std::sqrt(DistanceSquared(value));
}

template <typename T>
T Vector2<T>::DistanceSquared(const Vector2<T> &value) const noexcept
{
	T dx = x - value.x;
	T dy = y - value.y;
	return (dx * dx) + (dy * dy);
}

template <typename T>
void Vector2<T>::Normalize() noexcept
{
	T length = Length();
	if (length == T(0.0))
	{
		return;
	}

	T num = T(1.0) / length;
	x *= num;
	y *= num;
}

template <typename T>
bool Vector2<T>::IsZero() const noexcept
{
	return ((this->x == T(0.0)) && (this->y == T(0.0)));
}
template <typename T>
bool Vector2<T>::IsInfinity() const noexcept
{
	return ((this->x == std::numeric_limits<T>::max()) &&
		(this->y == std::numeric_limits<T>::max()));
}


template <typename T>
Vector2<T> Vector2<T>::Normalize(const Vector2<T> &vector) noexcept
{
	T length = vector.Length();
	if (length == T(0.0))
	{
		return vector;
	}

	T num = T(1.0) / length;	
	return Vector2<T>(vector.x * num, vector.y * num);
}

template <typename T>
Vector2<T> Vector2<T>::MakeOrtogonal(const Vector2<T> & v) noexcept
{
	return Vector2<T>(-v.y, v.x);
}



template <typename T>
T Vector2<T>::Distance(const Vector2<T> &value1, const Vector2<T> &value2) noexcept
{	
	return value1.Distance(value2);
}

template <typename T>
T Vector2<T>::DistanceSquared(const Vector2<T> &value1, const Vector2<T> &value2) noexcept
{
	return value1.DistanceSquared(value2);
}
	

template <typename T>
T Vector2<T>::CrossSquared(const Vector2<T> &left, const Vector2<T> &right) noexcept
{
	T aa = Vector2::Dot(left, left);
	T bb = Vector2::Dot(right, right);
	T ab = Vector2::Dot(left, right);
	return (aa * bb - ab * ab);
}

template <typename T>
Vector2<T> Vector2<T>::Center(const Vector2<T> &left, const Vector2<T> &right) noexcept
{	
	return Vector2((left.x + right.x) / 2.0f,
		(left.y + right.y) / 2.0f);
}

template <typename T>
T Vector2<T>::Dot(const Vector2<T> &left, const Vector2<T> &right) noexcept
{
	return (left.x * right.x + left.y * right.y);
}

template <typename T>
int Vector2<T>::SortCompare(const Vector2<T> & p1, const Vector2<T> & p2) noexcept
{
	if (p1.x > p2.x) return 1;
	if (p1.x < p2.x) return -1;

	// x coordintaes are the same, check y
	if (p1.y > p2.y) return 1;
	if (p1.y < p2.y) return -1;

	//y coordinates are the same, points are equal
	return 0;
};

#ifdef MYMATH_NAMESPACE 
#	define MY_MATH_NS MYMATH_NAMESPACE 
#else
#	define MY_MATH_NS 
#endif
	
template struct MY_MATH_NS::Vector2<float>;
template struct MY_MATH_NS::Vector2<double>;