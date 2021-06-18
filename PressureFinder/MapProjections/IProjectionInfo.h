#ifndef IPROJECTION_INFO_H
#define IPROJECTION_INFO_H

#include <functional>

#include "MapProjectionStructures.h"

//#define USE_VIRTUAL_INTERFACE

#ifdef USE_VIRTUAL_INTERFACE
#	define OVERRIDE override
#else
#	define OVERRIDE
#endif

namespace Projections 
{
	class IProjectionInfo
	{
	public:
		
		const PROJECTION curProjection;

		virtual ~IProjectionInfo() = default;

		virtual const char* GetName() const
		{
			return "";
		}

#ifdef USE_VIRTUAL_INTERFACE
		/*
		template <typename PixelType = int>
		virtual RET_VAL(PixelType, std::is_integral) Project(const Coordinate & c) const = 0;


		template <typename PixelType = float>
		virtual RET_VAL(PixelType, std::is_floating_point) Project(const Coordinate & c) const = 0;


		template <typename PixelType = int, bool Normalize = true>
		virtual Coordinate ProjectInverse(const Pixel<PixelType> & p) const = 0;

		template <typename InputProj>
		virtual void SetFrame(InputProj * proj, bool keepAR = true) = 0;

		template <typename InputProj>
		virtual void SetFrame(InputProj * proj, double w, double h, bool keepAR = true) = 0;
		*/

		virtual void SetFrame(const ProjectionFrame & frame) = 0;		
		virtual void SetFrame(const Coordinate & botLeft, const Coordinate & topRight, MyRealType w, MyRealType h, STEP_TYPE stepType, bool keepAR = true) = 0;
		virtual void SetRawFrame(const Coordinate & botLeft, const Coordinate & topRight, MyRealType w, MyRealType h, STEP_TYPE stepType, bool keepAR = true) = 0;
		virtual void SetFrameFromAABB(const Coordinate & min, const Coordinate & max, MyRealType w, MyRealType h, STEP_TYPE stepType, bool keepAR = true) = 0;
		
		
		virtual Coordinate GetTopLeftCorner() const = 0;
		virtual Coordinate GetDeltaStep() const = 0;
		virtual const ProjectionFrame & GetFrame() const = 0;
	
		virtual void LineBresenham(Pixel<int> start, Pixel<int> end,
			std::function<void(int x, int y)> callback) const = 0;

		virtual void ComputeAABB(Coordinate & min, Coordinate & max) const = 0;
#endif
	protected:
		IProjectionInfo(PROJECTION curProjection) : curProjection(curProjection) {};
	};

}

#endif
