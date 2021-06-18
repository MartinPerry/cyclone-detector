#ifndef PROJECTION_INFO_H
#define PROJECTION_INFO_H

#include <cstdint>
#include <cmath>
#include <vector>
#include <complex>
#include <functional>
#include <tuple>


#include <string>

#include "./GeoCoordinate.h"
#include "./MapProjectionStructures.h"
#include "./MapProjectionUtils.h"

#include "./IProjectionInfo.h"


#define RET_VAL(PixelType, enable_cond) \
	typename std::enable_if<enable_cond<PixelType>::value, Pixel<PixelType>>::type

namespace Projections
{
	template <typename Proj>
	class ProjectionInfo : public IProjectionInfo
	{
	public:				
		virtual ~ProjectionInfo() = default;

		const char* GetName() const;

		template <typename PixelType = int>
		RET_VAL(PixelType, std::is_integral) Project(const Coordinate & c) const;


		template <typename PixelType = float>
		RET_VAL(PixelType, std::is_floating_point) Project(const Coordinate & c) const;
        
		template <typename PixelType = int, bool Normalize = true>
		Coordinate ProjectInverse(const Pixel<PixelType> & p) const;

        
        
		template <typename InputProj>
		void SetFrame(InputProj * proj, bool keepAR = true);

		template <typename InputProj>
		void SetFrame(InputProj * proj, MyRealType w, MyRealType h, STEP_TYPE stepType, bool keepAR = true);

		void SetFrame(const ProjectionFrame & frame) OVERRIDE;		
		void SetFrame(const Coordinate & botLeft, const Coordinate & topRight, MyRealType w, MyRealType h, STEP_TYPE stepType, bool keepAR = true) OVERRIDE;
		void SetRawFrame(const Coordinate & botLeft, const Coordinate & topRight, MyRealType w, MyRealType h, STEP_TYPE stepType, bool keepAR = true) OVERRIDE;
		void SetFrameFromAABB(const Coordinate & min, const Coordinate & max, MyRealType w, MyRealType h, STEP_TYPE stepType, bool keepAR = true) OVERRIDE;


		Coordinate GetTopLeftCorner() const OVERRIDE;
		Coordinate GetDeltaStep() const OVERRIDE;
		const ProjectionFrame & GetFrame() const OVERRIDE;


		template <typename T = int>
		T GetFrameWidth() const { return static_cast<T>(this->frame.w); }
		template <typename T = int>
		T GetFrameHeight() const { return static_cast<T>(this->frame.h); }
		
		void LineBresenham(Pixel<int> start, Pixel<int> end, 
			std::function<void(int x, int y)> callback) const OVERRIDE;



		void ComputeAABB(Coordinate & min, Coordinate & max) const OVERRIDE;


	protected:

        ///<summary>
        /// Helper struct for actual projections return types
        ///</summary>
		struct ProjectedValueInverse
		{
			Latitude lat;
			Longitude lon;
		};
        
        ///<summary>
        /// Helper struct for actual projections return types
        ///</summary>
        struct ProjectedValue
        {
            MyRealType x;
            MyRealType y;
        };
		

		ProjectionFrame frame;

		std::tuple<double, double, double, double> 
			GetFrameBotLeftTopRight(const Coordinate & botLeft, const Coordinate & topRight);

		void CalculateWrapRepeat(const Coordinate& botLeft, const Coordinate& topRight);

		ProjectionInfo(PROJECTION curProjection);
	};



	/// <summary>
	/// Set current data active frame based on existing projection
	/// </summary>
	/// <param name="proj">existing projection</param>
	/// <param name="keepAR">keep AR of data (default: true) 
	/// if yes, data are enlarged and not 1:1 to bounding box to keep AR
	/// </param>
	template <typename Proj>
	template <typename InputProj>
	void ProjectionInfo<Proj>::SetFrame(InputProj * proj, bool keepAR)
	{
		auto f = proj->GetFrame();
		this->SetFrameFromAABB(f.min, f.max, f.w, f.h, f.stepType, keepAR);
	};


	/// <summary>
	/// Set current data active frame based on existing projection
	/// </summary>
	/// <param name="proj">existing projection</param>
	/// <param name="w">frame width</param>
	/// <param name="h">frame height</param>
	/// <param name="keepAR">keep AR of data (default: true) 
	/// if yes, data are enlarged and not 1:1 to bounding box to keep AR
	/// </param>
	template <typename Proj>
	template <typename InputProj>
	void ProjectionInfo<Proj>::SetFrame(InputProj * proj, MyRealType w, MyRealType h, STEP_TYPE stepType, bool keepAR)
	{
		auto f = proj->GetFrame();
		this->SetFrame(f.min, f.max, w, h, stepType, keepAR);
	};



	/// <summary>
	/// Project Coordinate point to pixel
	/// </summary>
	/// <param name="c"></param>
	/// <returns></returns>
	template <typename Proj>
	template <typename PixelType>	
	RET_VAL(PixelType, std::is_integral) ProjectionInfo<Proj>::Project(const Coordinate & c) const
	{

		ProjectedValue rawPixel = static_cast<const Proj*>(this)->ProjectInternal(c);

		Pixel<PixelType> p{
			static_cast<PixelType>(std::round(rawPixel.x * this->frame.wAR - this->frame.projPrecomX)),
			static_cast<PixelType>(std::round(-rawPixel.y * this->frame.hAR - this->frame.projPrecomY))
		};

		//rawPixel.x = rawPixel.x - this->frame.minPixelOffsetX;
		//rawPixel.y = rawPixel.y - this->frame.minPixelOffsetY;

		//Pixel<PixelType> p;
		//p.x = static_cast<PixelType>(std::round(this->frame.wPadding + (rawPixel.x * this->frame.wAR)));
		//p.y = static_cast<PixelType>(std::round(this->frame.h - this->frame.hPadding - (rawPixel.y * this->frame.hAR)));

		return p;
	};

	template <typename Proj>
	template <typename PixelType>	
	RET_VAL(PixelType, std::is_floating_point) ProjectionInfo<Proj>::Project(const Coordinate & c) const
	{

		//project value and get "pseudo" pixel coordinate
		ProjectedValue rawPixel = static_cast<const Proj*>(this)->ProjectInternal(c);

		Pixel<PixelType> p{
			static_cast<PixelType>(rawPixel.x * this->frame.wAR - this->frame.projPrecomX),
			static_cast<PixelType>(-rawPixel.y * this->frame.hAR - this->frame.projPrecomY)
		};

		//move our pseoude pixel to "origin"
		//rawPixel.x = rawPixel.x - this->frame.minPixelOffsetX;
		//rawPixel.y = rawPixel.y - this->frame.minPixelOffsetY;

		//calculate pixel in final frame
		//Pixel<PixelType> p;
		//p.x = static_cast<PixelType>(this->frame.wPadding + (rawPixel.x * this->frame.wAR));
		//p.y = static_cast<PixelType>(this->frame.h - this->frame.hPadding - (rawPixel.y * this->frame.hAR));

		return p;
	};

	/// <summary>
	/// Project pixel to coordinate
	/// </summary>
	/// <param name="p"></param>
	/// <returns></returns>
	template <typename Proj>
	template <typename PixelType, bool Normalize>
	Coordinate ProjectionInfo<Proj>::ProjectInverse(const Pixel<PixelType> & p) const
	{

		//double xx = (static_cast<double>(p.x) - this->frame.wPadding + this->frame.wAR * this->frame.minPixelOffset.x);
		MyRealType xx = (static_cast<MyRealType>(p.x) + this->frame.projPrecomX);
		xx /= this->frame.wAR;

		//double yy = (static_cast<double>(p.y) - this->frame.h + this->frame.hPadding - this->frame.hAR * this->frame.minPixelOffset.y);
		MyRealType yy = (static_cast<MyRealType>(p.y) + this->frame.projPrecomY);
		yy /= -this->frame.hAR;


		
		ProjectedValueInverse pi = static_cast<const Proj*>(this)->ProjectInverseInternal(xx, yy);

		Coordinate c;
		c.lat = pi.lat;
		c.lon = pi.lon;

		if (Normalize)
		{
			c.lat.Normalize();
			c.lon.Normalize();			
		}
		
		return c;
	};

	
	

};

#endif
