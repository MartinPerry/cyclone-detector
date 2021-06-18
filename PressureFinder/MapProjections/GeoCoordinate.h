#ifndef GEOCORDINATE_H
#define GEOCORDINATE_H

#include <cmath>

typedef double MyRealType;

struct AngleUtils
{
	static MyRealType radToDeg(MyRealType val) { return val * MyRealType(57.2957795); }
	static MyRealType degToRad(MyRealType val) { return val * MyRealType(0.0174532925); }
};

template <typename T>
struct IAngle
{
	IAngle() : valRad(0), valDeg(0) {};
	static T deg(MyRealType val) { return T(AngleUtils::degToRad(val), val); };
	static T rad(MyRealType val) { return T(val, AngleUtils::radToDeg(val)); };
	
	inline MyRealType deg() const { return valDeg; };
	inline MyRealType rad() const { return valRad; };

	inline T operator -() { return T(-valRad, -valDeg); };

protected:
	IAngle(MyRealType valRad, MyRealType valDeg) : valRad(valRad), valDeg(valDeg) {};
	MyRealType valRad;
	MyRealType valDeg;
};

struct AngleValue : public IAngle<AngleValue>
{
	AngleValue() : IAngle() {};

	friend struct IAngle<AngleValue>;

protected:
	AngleValue(MyRealType valRad, MyRealType valDeg) : IAngle(valRad, valDeg) {};
};

struct Latitude : public IAngle<Latitude>
{
	Latitude() : IAngle() {};
	Latitude(const Latitude & a) : IAngle(a.rad(), a.deg()) {};
	Latitude(const AngleValue & a) : IAngle(a.rad(), a.deg()) {};

	/// <summary>
	/// Normlization will only clamp latitude to [-90, 90] deg interval
	/// since Latitude has no wrap around
	/// </summary>
	void Normalize()
	{
		this->Clamp();
		//valDeg = (valDeg > 90) ? (valDeg - 180) : valDeg;		
		valRad = Latitude::deg(valDeg).rad();
	};

	/// <summary>
	/// Clamp to [-90, 90] interval
	/// </summary>
	void Clamp()
	{
		//valDeg = (valDeg > 90) ? (valDeg - 180) : valDeg;
		valDeg = (valDeg > 90) ? 90.0 : ((valDeg < -90) ? -90.0 : valDeg);
		valRad = Latitude::deg(valDeg).rad();
	};

	friend struct IAngle<Latitude>;
protected:
	Latitude(MyRealType valRad, MyRealType valDeg) : IAngle(valRad, valDeg) {};
};

struct Longitude : public IAngle<Longitude>
{
	Longitude() : IAngle() {};
	Longitude(const Longitude & a) : IAngle(a.rad(), a.deg()) {};
	Longitude(const AngleValue & a) : IAngle(a.rad(), a.deg()) {};

	/// <summary>
	/// Normlize to [-180, 180] interval with a wrap-around
	/// </summary>
	void Normalize()
	{
		while (valDeg < -180)
		{
			valDeg += 360;
		}
		while (valDeg > 180)
		{
			valDeg -= 360;
		}
		//valDeg = std::fmod(std::fmod(valDeg, 360) + 540, 360) - 180;		
		valRad = Longitude::deg(valDeg).rad();
	};

	/// <summary>
	/// Clamp to [-180, 180] interval
	/// </summary>
	void Clamp()
	{
		valDeg = (valDeg > 180.0) ? 180.0 : ((valDeg < -180.0) ? -180.0 : valDeg);
		valRad = Latitude::deg(valDeg).rad();
	};

	friend struct IAngle<Longitude>;
protected:
	Longitude(MyRealType valRad, MyRealType valDeg) : IAngle(valRad, valDeg) {};
};


//==============================================================
//String literall operator for Angle only
//Latitude and longitude can be created from Angle

inline AngleValue operator "" _deg(long double value)
{
	return AngleValue::deg(value);
}

inline AngleValue operator "" _rad(long double value)
{
	return  AngleValue::rad(value);
}




#endif
