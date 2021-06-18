#ifndef REPROJECTION_H
#define REPROJECTION_H

#include <vector>


#include "./MapProjectionStructures.h"
#include "./ProjectionInfo.h"

namespace Projections
{

	/// <summary>
	/// Reprojection structure
	/// It holds info needed to reproject data from one projection
	/// to another
	/// pixels - reprojection info
	/// pixels[to] = from
	/// 
	/// Template type is type of reprojection Pixel
	/// By default its int -> reprojection can be in range of int
	/// However in many cases we dont have such big images
	/// and we can use short
	/// 
	/// Calculates mapping: toData[index] = fromData[reprojection[index]]
	/// </summary>
	/*
		template <typename T = int,
			typename = typename std::enable_if<
			std::is_same<T, int>::value ||
			std::is_same<T, short>::value>::type
		>
	*/
	template <typename T = int>
	struct Reprojection
	{	
		int inW;
		int inH;
		int outW;
		int outH;
		std::vector<Pixel<T>> pixels; //[to] = from

		Reprojection() : 
			inW(0),
			inH(0),
			outW(0),
			outH(0)
		{
		}

		/// <summary>
		/// Create cache name in format:
		/// reproj_from_w_h_to_w_h
		/// </summary>
		/// <param name="from"></param>
		/// <param name="to"></param>
		/// <returns></returns>
		template <typename FromProjection, typename ToProjection>
		static std::string CreateCacheName(FromProjection* from, ToProjection* to)
		{
			std::string tmp = "reproj_";
			tmp += from->GetName();
			tmp += "_";
#ifdef USE_VIRTUAL_INTERFACE
			tmp += std::to_string(static_cast<int>(from->GetFrame().w));
			tmp += "_";
			tmp += std::to_string(static_cast<int>(from->GetFrame().h));
			tmp += "_";
#endif
			tmp += to->GetName();
#ifdef USE_VIRTUAL_INTERFACE
			tmp += "_";
			tmp += std::to_string(static_cast<int>(to->GetFrame().w));
			tmp += "_";
			tmp += std::to_string(static_cast<int>(to->GetFrame().h));
#endif
			
			return tmp;
		}

		/// <summary>
		/// Load reprojection from file		
		/// </summary>
		/// <param name="imProj"></param>
		/// <returns></returns>	
		static Reprojection<T> CreateFromFile(const std::string& fileName);
		
		/// <summary>
		/// Re-project data from -> to
		/// Calculates mapping: toData[index] = fromData[reprojection[index]]
		/// </summary>
		/// <param name="imProj"></param>
		/// <returns></returns>	
		template <typename FromProjection, typename ToProjection>
		static Reprojection<T> CreateReprojection(FromProjection* from, ToProjection* to)
		{

			Reprojection<T> reprojection;
			reprojection.pixels.resize(to->GetFrameWidth() * to->GetFrameHeight(), { -1, -1 });

			const auto& f = to->GetFrame();

			if ((f.repeatNegCount != 0) || (f.repeatPosCount != 0))
			{
				//we have multiple wrap around of the world

				//calculate full size of from projection image
				Projections::Coordinate bbMin, bbMax;

				bbMin.lat = -90.0_deg;
				bbMin.lon = -180.0_deg;

				bbMax.lat = 90.06_deg;
				bbMax.lon = 180.0_deg;

				Pixel<MyRealType> pp1 = from->template Project<MyRealType>(bbMin);
				Pixel<MyRealType> pp2 = from->template Project<MyRealType>(bbMax);

				MyRealType ww = pp2.x - pp1.x;
				MyRealType hh = pp1.y - pp2.y;
				
				for (int y = 0; y < to->GetFrameHeight(); y++)
				{
					int yw = y * to->GetFrameWidth();
					for (int x = 0; x < to->GetFrameWidth(); x++)
					{
						Pixel<T> p = Reprojection<T>::ReProject<int, T>({ x, y }, from, to);

						if ((p.x >= 0) &&
							(p.y >= 0) &&
							(p.x < from->GetFrameWidth()) &&
							(p.y < from->GetFrameHeight()))
						{
							reprojection.pixels[x + yw] = p;
						}


						int offset = static_cast<int>(ww);

						int px = p.x;

						MyRealType nc = f.repeatNegCount;						
						while (nc > 0)
						{
							p.x += offset;

							if ((p.x >= 0) &&
								(p.y >= 0) &&
								(p.x < from->GetFrameWidth()) &&
								(p.y < from->GetFrameHeight()))
							{
								reprojection.pixels[x + yw] = p;
							}
							nc--;
						}

						p.x = px; //restore p.x

						MyRealType pc = f.repeatPosCount;						
						while (pc > 0)
						{
							p.x -= offset;

							if ((p.x >= 0) &&
								(p.y >= 0) &&
								(p.x < from->GetFrameWidth()) &&
								(p.y < from->GetFrameHeight()))
							{
								reprojection.pixels[x + yw] = p;
							}
							pc--;
						}

					}
				}

			}						
			else if ((from->INDEPENDENT_LAT_LON) && (to->INDEPENDENT_LAT_LON))
			{
				//if x and y are independent, simplify

				std::vector<T> cacheX;
				cacheX.resize(to->GetFrameWidth());

				std::vector<T> cacheY;
				cacheY.resize(to->GetFrameHeight());

				for (int x = 0; x < to->GetFrameWidth(); x++)
				{
					Projections::Pixel<T> p = Reprojection<T>::ReProject<int, T>({ x, 0 }, from, to);
					cacheX[x] = p.x;
				}

				for (int y = 0; y < to->GetFrameHeight(); y++)
				{
					Projections::Pixel<T> p = Reprojection<T>::ReProject<int, T>({ 0, y }, from, to);
					cacheY[y] = p.y;

				}

				for (int y = 0; y < to->GetFrameHeight(); y++)
				{
					int yw = y * to->GetFrameWidth();
					for (int x = 0; x < to->GetFrameWidth(); x++)
					{
						Pixel<T> p;
						p.x = cacheX[x];
						p.y = cacheY[y];

						if (p.x < 0) continue;
						if (p.y < 0) continue;
						if (p.x >= from->GetFrameWidth()) continue;
						if (p.y >= from->GetFrameHeight()) continue;

						reprojection.pixels[x + yw] = p;

					}
				}
			}			
			else 
			{				
				for (int y = 0; y < to->GetFrameHeight(); y++)
				{
					int yw = y * to->GetFrameWidth();
					for (int x = 0; x < to->GetFrameWidth(); x++)
					{
						Pixel<T> p = Reprojection<T>::ReProject<int, T>({ x, y }, from, to);

						if ((p.x >= 0) && 
							(p.y >= 0) && 
							(p.x < from->GetFrameWidth()) &&
							(p.y < from->GetFrameHeight()))
						{
							reprojection.pixels[x + yw] = p;
						}						
					}
				}
			}

			reprojection.inW = from->GetFrameWidth();
			reprojection.inH = from->GetFrameHeight();
			reprojection.outW = to->GetFrameWidth();
			reprojection.outH = to->GetFrameHeight();

			return reprojection;
		};

		
		
		/// <summary>
		/// Save reprojection to file
		/// </summary>
		/// <param name="fileName"></param>
		void SaveToFile(const std::string& fileName);


		/// <summary>
		/// Reproject inputData based on reproj with Nerest Neighbor interpolation.
		/// Output array has size reproj.outW * reproj.outH
		/// Output array must be released with delete[]
		/// 
		/// Template parameters:
		/// DataType - type of input data		
		/// Out - output structure - can be raw array of std::vector
		/// ChannelsCount - number of channels in input / output data
		/// </summary>
		/// <param name="reproj"></param>
		/// <param name="inputData"></param>
		/// <param name="NO_VALUE"></param>
		/// <returns></returns>
		template <typename DataType, typename Out = DataType*, size_t ChannelsCount = 1>
		Out ReprojectDataNerestNeighbor(DataType* inputData, const DataType NO_VALUE) const
		{
			size_t count = this->outW * this->outH;

			Out output;

			if constexpr (std::is_same<Out, DataType*>::value)
			{
				output = new DataType[count * ChannelsCount];
			}
			else if constexpr (std::is_same<Out, std::vector<DataType>>::value)
			{
				output.resize(count * ChannelsCount);
			}

			for (size_t index = 0; index < count; index++)
			{
				T x = static_cast<int>(this->pixels[index].x);
				T y = static_cast<int>(this->pixels[index].y);

				if ((x == -1) || (y == -1))
				{
					//outside of the model - no data - put there NO_VALUE
					if constexpr (ChannelsCount == 1)
					{
						output[index] = NO_VALUE;
					}
					else
					{
						for (size_t i = 0; i < ChannelsCount; i++)
						{
							output[index * ChannelsCount + i] = NO_VALUE;
						}
					}
				}
				else
				{
					size_t origIndex = x + y * this->inW;
					if constexpr (ChannelsCount == 1)
					{
						output[index] = inputData[origIndex];
					}
					else
					{
						for (size_t i = 0; i < ChannelsCount; i++)
						{
							output[index * ChannelsCount + i] = inputData[origIndex * ChannelsCount + i];
						}
					}
				}
			}

			return output;
		}

		/// <summary>
		/// Reproject inputData based on reproj with Bilinear interpolation.
		/// Note: Usable only if T is nor int number
		/// 
		/// Output array has size reproj.outW * reproj.outH
		/// Output array must be released with delete[]
		/// 
		/// Template parameters:
		/// DataType - type of input data		
		/// Out - output structure - can be raw array of std::vector
		/// ChannelsCount - number of channels in input / output data
		/// </summary>
		/// <param name="reproj"></param>
		/// <param name="inputData"></param>
		/// <param name="NO_VALUE"></param>
		/// <returns></returns>
		template <typename DataType, typename Out = DataType*, size_t ChannelsCount = 1>
		Out ReprojectDataBilinear(DataType* inputData, const DataType NO_VALUE) const
		{
			size_t count = this->outW * this->outH;

			Out output;

			if constexpr (std::is_same<Out, DataType*>::value)
			{
				output = new DataType[count * ChannelsCount];
			}
			else if constexpr (std::is_same<Out, std::vector<DataType>>::value)
			{
				output.resize(count * ChannelsCount);
			}
			
			for (size_t index = 0; index < count; index++)
			{
				T x = this->pixels[index].x;
				T y = this->pixels[index].y;

				if ((x == -1) || (y == -1))
				{
					//outside of the model - no data - put there NO_VALUE
					if constexpr (ChannelsCount == 1)
					{
						output[index] = NO_VALUE;
					}
					else
					{
						for (size_t i = 0; i < ChannelsCount; i++)
						{
							output[index * ChannelsCount + i] = NO_VALUE;
						}
					}
				}
				else
				{		
					//no floor, just cast, because values x and y are non-negative
					int px = static_cast<int>(x);
					int py = static_cast<int>(y);

					double tx = x - px;
					double ty = y - py;

					int x1p = (px + 1 >= this->inW) ? this->inW - 1 : px + 1;
					int y1p = (py + 1 >= this->inH) ? this->inH - 1 : py + 1;

					const DataType* c00 = &inputData[(px + py * this->inW) * ChannelsCount];
					const DataType* c10 = &inputData[(x1p + py * this->inW) * ChannelsCount];
					const DataType* c01 = &inputData[(px + y1p * this->inW) * ChannelsCount];
					const DataType* c11 = &inputData[(x1p + y1p * this->inW) * ChannelsCount];



					for (size_t i = 0; i < ChannelsCount; i++)
					{
						auto a = c00[i] * (1 - tx) + c10[i] * tx;
						auto b = c01[i] * (1 - tx) + c11[i] * tx;
						auto res = a * (1 - ty) + b * ty;
						
						output[index * ChannelsCount + i] = res;						
					}
				}
			}


			return output;
		}


		/// <summary>
		/// Reproject inputData based on reproj with Bicubic interpolation.
		/// Note: Usable only if T is nor int number
		/// 
		/// Output array has size reproj.outW * reproj.outH
		/// Output array must be released with delete[]
		/// 
		/// Template parameters:
		/// DataType - type of input data		
		/// Out - output structure - can be raw array of std::vector
		/// ChannelsCount - number of channels in input / output data
		/// </summary>
		/// <param name="reproj"></param>
		/// <param name="inputData"></param>
		/// <param name="NO_VALUE"></param>
		/// <returns></returns>
		template <typename DataType, typename Out = DataType*, size_t ChannelsCount = 1>
		Out ReprojectDataBicubic(DataType* inputData, const DataType NO_VALUE) const
		{
			size_t count = this->outW * this->outH;

			Out output;

			if constexpr (std::is_same<Out, DataType*>::value)
			{
				output = new DataType[count * ChannelsCount];
			}
			else if constexpr (std::is_same<Out, std::vector<DataType>>::value)
			{
				output.resize(count * ChannelsCount);
			}

			for (size_t index = 0; index < count; index++)
			{
				T x = this->pixels[index].x;
				T y = this->pixels[index].y;

				if ((x == -1) || (y == -1))
				{
					//outside of the model - no data - put there NO_VALUE
					if constexpr (ChannelsCount == 1)
					{
						output[index] = NO_VALUE;
					}
					else
					{
						for (size_t i = 0; i < ChannelsCount; i++)
						{
							output[index * ChannelsCount + i] = NO_VALUE;
						}
					}
				}
				else
				{
					//no floor, just cast, because values x and y are non-negative
					int px = static_cast<int>(x);
					int py = static_cast<int>(y);

					double fx = x - px;
					double fy = y - py;
				

					//we'll need the second and third powers
					//of f to compute our filter weights
					double f2x = fx * fx;
					double f3x = f2x * fx;

					double f2y = fy * fy;
					double f3y = f2y * fy;
					
					double fmpF1x = (1.0f - fx);
					double f12x = fmpF1x * fmpF1x;
					double f13x = f12x * fmpF1x;

					double fmpF1y = (1.0f - fy);
					double f12y = fmpF1y * fmpF1y;
					double f13y = f12y * fmpF1y;


					//compute the filter weights

					double w0x = (f13x);
					double w1x = (4.0f + 3.0f * f3x - 6.0f * f2x);
					double w2x = (4.0f + 3.0f * f13x - 6.0f * f12x);
					double w3x = (f3x);

					double w0y = (f13y);
					double w1y = (4.0f + 3.0f * f3y - 6.0f * f2y);
					double w2y = (4.0f + 3.0f * f13y - 6.0f * f12y);
					double w3y = (f3y);


					int x1m = (px < 1) ? 0 : px - 1;
					int x1p = (px + 1 >= this->inW) ? this->inW - 1 : px + 1;
					int x2p = (px + 2 >= this->inW) ? this->inW - 2 : px + 2;

					int y1m = (py < 1) ? 0 : py - 1;
					int y1p = (py + 1 >= this->inH) ? this->inH - 1 : py + 1;
					int y2p = (py + 2 >= this->inH) ? this->inH - 2 : py + 2;

					
					const DataType* p00 = &inputData[(x1m + y1m * this->inW) * ChannelsCount];
					const DataType* p10 = &inputData[(px + y1m * this->inW) * ChannelsCount];
					const DataType* p20 = &inputData[(x1p + y1m * this->inW) * ChannelsCount];
					const DataType* p30 = &inputData[(x2p + y1m * this->inW) * ChannelsCount];

					const DataType* p01 = &inputData[(x1m + py * this->inW) * ChannelsCount];
					const DataType* p11 = &inputData[(px + py * this->inW) * ChannelsCount];
					const DataType* p21 = &inputData[(x1p + py * this->inW) * ChannelsCount];
					const DataType* p31 = &inputData[(x2p + py * this->inW) * ChannelsCount];

					const DataType* p02 = &inputData[(x1m + y1p * this->inW) * ChannelsCount];
					const DataType* p12 = &inputData[(px + y1p * this->inW) * ChannelsCount];
					const DataType* p22 = &inputData[(x1p + y1p * this->inW) * ChannelsCount];
					const DataType* p32 = &inputData[(x2p + y1p * this->inW) * ChannelsCount];

					const DataType* p03 = &inputData[(x1m + y2p * this->inW) * ChannelsCount];
					const DataType* p13 = &inputData[(px + y2p * this->inW) * ChannelsCount];
					const DataType* p23 = &inputData[(x1p + y2p * this->inW) * ChannelsCount];
					const DataType* p33 = &inputData[(x2p + y2p * this->inW) * ChannelsCount];


					for (size_t i = 0; i < ChannelsCount; i++)
					{
												
						double res = (1.0 / 36.0) * (
							w0y * (p00[i] * w0x
								+ p10[i] * w1x
								+ p20[i] * w2x
								+ p30[i] * w3x)

							+ w1y * (p01[i] * w0x
								+ p11[i] * w1x
								+ p21[i] * w2x
								+ p31[i] * w3x)

							+ w2y * (p02[i] * w0x
								+ p12[i] * w1x
								+ p22[i] * w2x
								+ p32[i] * w3x)

							+ w3y * (p03[i] * w0x
								+ p13[i] * w1x
								+ p23[i] * w2x
								+ p33[i] * w3x)
							);

						output[index * ChannelsCount + i] = static_cast<DataType>(res);
					}
				}
			}


			return output;
		}

		/// <summary>
		/// Reproject single pixel from -> to
		/// </summary>
		/// <param name="p"></param>
		/// <param name="from"></param>
		/// <param name="to"></param>
		/// <returns></returns>
		template <typename InPixelType, typename OutPixelType,
			typename FromProjection, typename ToProjection>
			static Pixel<OutPixelType> ReProject(Pixel<InPixelType> p,
				const FromProjection* from, const ToProjection* to)
		{
			Coordinate cc = to->template ProjectInverse<InPixelType, false>(p);
			return from->template Project<OutPixelType>(cc);
		};
	};

}

#endif
