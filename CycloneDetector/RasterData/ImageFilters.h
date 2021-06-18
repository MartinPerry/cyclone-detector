#ifndef IMAGE_FILTERS_H
#define IMAGE_FILTERS_H

template <typename T>
class Image2d;

template <typename T>
class ImageGausFilter;

struct ConvolutionFilter2d;
struct ConvolutionFilterSeparable;

#include <vector>
#include <array>
#include <functional>

#ifdef ENABLE_SIMD
#	include <immintrin.h>
#endif

#include "./ImageUtils.h"

//=================================================================================================

#ifndef CALL_MEMBER_FN
#	define CALL_MEMBER_FN(object, ptrToMember)  ((object).*(ptrToMember))
#endif

struct ConvolutionFilter
{
	template <typename T, typename F>
	using RunCallbackPtr = float(F::*)(const Image2d<T> &, int, int) const;


	static void ResultUpdateNone(float & res);
	static void ResultUpdateSqr(float & res);

	NeighborhoodKernel n;

	ImageUtils::BorderMode borderMode;

	void(*filterResultCallback)(float & res); //callback to update result value after filter finish run
											  //default: ResultUpdateNone		

	virtual ConvolutionFilter2d * GetAsConvolutionFilter2d();
	virtual ConvolutionFilterSeparable * GetAsConvolutionFilterSeparable();

	virtual const ConvolutionFilter2d * GetAsConvolutionFilter2d() const;
	virtual const ConvolutionFilterSeparable * GetAsConvolutionFilterSeparable() const;

	int GetSize() const;

	void SetActiveChannel(size_t channel) const;

	template <typename T, typename V, typename F>
	static void LoopHelper(const Image2d<T> & input, Image2d<V> & output,
		int yStart, int yEnd,
		int xStart, int xEnd,
		const F & f, ConvolutionFilter::RunCallbackPtr<T, F> p);

	//smazat
	static bool enableSimd;

protected:
	ConvolutionFilter(NeighborhoodKernel k, ImageUtils::BorderMode borderMode,
		void(*filterResultCallback)(float & res));

	mutable size_t activeChannel;
};

//=================================================================================================

struct ConvolutionFilterSeparable : public ConvolutionFilter
{
	enum class Order
	{
		X_Y,
		Y_X
	};

	static ConvolutionFilterSeparable CreateSobelX();
	static ConvolutionFilterSeparable CreateSobelY();
	static ConvolutionFilterSeparable CreateScharrX();
	static ConvolutionFilterSeparable CreateScharrY();

	float mulX;
	float mulY;

	std::vector<float> kernelDataX;
	std::vector<float> kernelDataY;

	Order order;

	ConvolutionFilterSeparable();
	ConvolutionFilterSeparable(const std::vector<float> & data);
	ConvolutionFilterSeparable(const std::vector<float> & dataX,
		const std::vector<float> & dataY, Order order);

	ConvolutionFilterSeparable * GetAsConvolutionFilterSeparable() override;
	const ConvolutionFilterSeparable * GetAsConvolutionFilterSeparable() const override;

	template <typename T>
	float RunX(const Image2d<T> & input, int x, int y) const;

	template <typename T>
	float RunY(const Image2d<T> & input, int x, int y) const;

	template <typename T>
	float RunWithBorderX(const Image2d<T> & input, int x, int y) const;

	template <typename T>
	float RunWithBorderY(const Image2d<T> & input, int x, int y) const;


	template <typename T, typename V>
	void LoopHelperSingleChannelX(const Image2d<T> & input, Image2d<V> & output,
		int yStart, int yEnd,
		int xStart, int xEnd) const;

	template <typename T, typename V>
	void LoopHelperSingleChannelY(const Image2d<T> & input, Image2d<V> & output,
		int yStart, int yEnd,
		int xStart, int xEnd) const;
};

//=================================================================================================

struct ConvolutionFilter2d : public ConvolutionFilter
{
	static ConvolutionFilter2d CreatePrevitX();
	static ConvolutionFilter2d CreatePrevitY();
	static ConvolutionFilter2d CreateSobelX();
	static ConvolutionFilter2d CreateSobelY();
	static ConvolutionFilter2d CreateScharrX();
	static ConvolutionFilter2d CreateScharrY();
	static ConvolutionFilter2d CreateSharpen();

	static ConvolutionFilter2d CreateDxx();
	static ConvolutionFilter2d CreateDyy();
	static ConvolutionFilter2d CreateDxy();

	static ConvolutionFilter2d CreateSquare(int radius, float val);
	static ConvolutionFilter2d CreateCircle(int radius, float val);

	float mul;
	std::vector<float> kernelData;

	static ConvolutionFilter2d CreateReflected(const ConvolutionFilter2d & f);

	ConvolutionFilter2d();
	ConvolutionFilter2d(const std::vector<float> & data);

	ConvolutionFilter2d * GetAsConvolutionFilter2d() override;
	const ConvolutionFilter2d * GetAsConvolutionFilter2d() const override;

	float * GetValueAt(int x, int y);
	const float * GetValueAt(int x, int y) const;

	void Negate();
	void Multiply(float val);

	template <typename T>
	float Run(const Image2d<T> & input, int x, int y) const;

	template <typename T>
	float RunWithBorder(const Image2d<T> & input, int x, int y) const;

protected:

};

//=================================================================================================

template <typename T>
class ImageFilters
{
public:

	static void ThressholdBinary(const Image2d<T> & input, Image2d<T> & output,
		const T * thressholds, const std::array<T, 2> * binaryValues);

	static void ThressholdAdaptive(const Image2d<T> & input, Image2d<T> & output,
		ImageGausFilter<T> & gauss, const std::array<T, 2> * binaryValues, double delta);

	static void EdgeDetection(const Image2d<T> & input, Image2d<T> & output);

	static void EdgeDetection(const Image2d<T> & input, Image2d<T> & output,
		std::function<T(float, float)> edgeCallback);

	static void MedianBinaryImage(const Image2d<T> & input, Image2d<T> & output,
		NeighborhoodKernel k,
		ImageUtils::BorderMode borderMode = ImageUtils::BorderMode::CLAMP);

	template <typename V = T>
	static void RunConvolutionFilter(const Image2d<T> & input, Image2d<V> & output,
		const ConvolutionFilter2d & f);

	template <typename V = T>
	static void RunConvolutionFilter(const Image2d<T> & input, Image2d<V> & output,
		const ConvolutionFilterSeparable & f);

private:
	//temporary image used in Separable calculation
	//Image2d<T> temp;

};


/// <summary>
/// Helper method for convolution filter iteration to "safe" space
/// by writing the for-loops
/// </summary>
/// <param name="input"></param>
/// <param name="output"></param>
/// <param name="yStart"></param>
/// <param name="yEnd"></param>
/// <param name="xStart"></param>
/// <param name="xEnd"></param>
/// <param name="f"></param>
/// <param name="p"></param>
template <typename T, typename V, typename F>
void ConvolutionFilter::LoopHelper(const Image2d<T> & input, Image2d<V> & output,
	int yStart, int yEnd,
	int xStart, int xEnd,
	const F & f, ConvolutionFilter::RunCallbackPtr<T, F> p)
{
	for (int y = yStart; y < yEnd; y++)
	{
		size_t yW = size_t(y) * output.GetWidth();
		for (int x = xStart; x < xEnd; x++)
		{
			float tmp = CALL_MEMBER_FN(f, p)(input, x, y);
			output.SetValue(tmp, f.activeChannel, size_t(x) + yW);
		}
	}
};



template <typename T, typename V>
void ConvolutionFilterSeparable::LoopHelperSingleChannelX(const Image2d<T> & input, Image2d<V> & output,
	int yStart, int yEnd,
	int xStart, int xEnd) const
{
	if (this->kernelDataX.size() < 4)
	{
		//kernel is too small to be vectorized
		ConvolutionFilter::LoopHelper(input, output,
			yStart, yEnd,
			xStart, xEnd,
			*this, &ConvolutionFilterSeparable::RunX);
		return;
	}

	V * outData = output.GetData().data();

#ifdef ENABLE_SIMD
	const T * inData = input.GetData().data();
	
	int kn = this->n.GetSize();
	int kn4 = kn - kn % 4;
#endif

	for (int y = yStart; y < yEnd; y++)
	{
		size_t yW = size_t(y) * output.GetWidth();
		size_t outX = size_t(xStart) + yW;

		int xEndSimd = xStart;
#ifdef ENABLE_SIMD

		int len = xEnd - xStart;
		len = len - len % 4;
		xEndSimd = xStart + len;

		for (int x = xStart; x < xEndSimd; x += 4)
		{
			size_t inYW = size_t(y) * input.GetWidth();
			size_t x0 = ((x + 0) - this->n.radius) + inYW;
			size_t x1 = x0 + 1;
			size_t x2 = x0 + 2;
			size_t x3 = x0 + 3;

			__m128 sumDot = _mm_setzero_ps();

			int i = 0;
			for (; i < kn4; i += 4)
			{
				__m128 kx = _mm_set_ps1(this->kernelDataX[i + 0]);
				__m128 ky = _mm_set_ps1(this->kernelDataX[i + 1]);
				__m128 kz = _mm_set_ps1(this->kernelDataX[i + 2]);
				__m128 kw = _mm_set_ps1(this->kernelDataX[i + 3]);


				__m128 dx, dy, dz, dw;

				if constexpr (std::is_same<T, uint8_t>::value)
				{
					//we need co convert uint8_t inputs to float
					__m128i u8_0 = _mm_loadu_si128((const __m128i*)(inData + x0));
					__m128i u8_1 = _mm_loadu_si128((const __m128i*)(inData + x1));
					__m128i u8_2 = _mm_loadu_si128((const __m128i*)(inData + x2));
					__m128i u8_3 = _mm_loadu_si128((const __m128i*)(inData + x3));

					__m128i u32_0 = _mm_unpacklo_epi16(
						_mm_unpacklo_epi8(u8_0, _mm_setzero_si128()),
						_mm_setzero_si128());
					__m128i u32_1 = _mm_unpacklo_epi16(
						_mm_unpacklo_epi8(u8_1, _mm_setzero_si128()),
						_mm_setzero_si128());
					__m128i u32_2 = _mm_unpacklo_epi16(
						_mm_unpacklo_epi8(u8_2, _mm_setzero_si128()),
						_mm_setzero_si128());
					__m128i u32_3 = _mm_unpacklo_epi16(
						_mm_unpacklo_epi8(u8_3, _mm_setzero_si128()),
						_mm_setzero_si128());

					dx = _mm_cvtepi32_ps(u32_0);
					dy = _mm_cvtepi32_ps(u32_1);
					dz = _mm_cvtepi32_ps(u32_2);
					dw = _mm_cvtepi32_ps(u32_3);

				}
				else
				{
					/*
					//load 8 consecutive values
					auto dd = _mm256_loadu_ps(inData + x0);

					//extract parts by shifting and casting to 4 values float
					dx = _mm256_castps256_ps128(dd);
					dy = _mm256_castps256_ps128(_mm256_permutevar8x32_ps(dd, _mm256_set_epi32(0, 0, 0, 0, 4, 3, 2, 1)));
					dz = _mm256_castps256_ps128(_mm256_permutevar8x32_ps(dd, _mm256_set_epi32(0, 0, 0, 0, 5, 4, 3, 2)));
					dw = _mm256_castps256_ps128(_mm256_permutevar8x32_ps(dd, _mm256_set_epi32(0, 0, 0, 0, 6, 5, 4, 3)));
					*/

					dx = _mm_loadu_ps(inData + x0);
					dy = _mm_loadu_ps(inData + x1);
					dz = _mm_loadu_ps(inData + x2);
					dw = _mm_loadu_ps(inData + x3);
				}

				//calculate 4 dots at once
				//[dx, dy, dz, dw] <dot> [kx, ky, kz, kw]

				auto mx = _mm_mul_ps(dx, kx); //dx * kx
				auto my = _mm_fmadd_ps(dy, ky, mx); //mx + dy * ky
				auto mz = _mm_fmadd_ps(dz, kz, my); //my + dz * kz
				auto res = _mm_fmadd_ps(dw, kw, mz); //mz + dw * kw

				sumDot = _mm_add_ps(sumDot, res);

				x0 += 4;
				x1 += 4;
				x2 += 4;
				x3 += 4;
			}

			for (; i < kn; i++)
			{
				auto v = _mm_set_ps1(this->kernelDataX[i]);
				auto v2 = _mm_set_ps(
					*(inData + x3), *(inData + x2),
					*(inData + x1), *(inData + x0)
				);

				sumDot = _mm_add_ps(sumDot, _mm_mul_ps(v, v2));

				x0++;
				x1++;
				x2++;
				x3++;
			}

			sumDot = _mm_mul_ps(sumDot, _mm_set_ps1(this->mulX));

			if constexpr (std::is_same<V, uint8_t>::value)
			{
				__m128i asInt = _mm_cvtps_epi32(sumDot);

				asInt = _mm_packus_epi32(asInt, asInt);	// Pack down to 16 bits
				asInt = _mm_packus_epi16(asInt, asInt); // Pack down to 8 bits

				uint32_t res = _mm_cvtsi128_si32(asInt); // Store the lower 32 bits	

				((uint32_t *)(outData + outX))[0] = res;
				outX += 4;
			}
			else
			{
				float tmpRes[4];
				_mm_store_ps(tmpRes, sumDot);

				outData[outX + 0] = tmpRes[0];
				outData[outX + 1] = tmpRes[1];
				outData[outX + 2] = tmpRes[2];
				outData[outX + 3] = tmpRes[3];
				outX += 4;
			}

		}
#endif		

		for (int x = xEndSimd; x < xEnd; x++)
		{
			float tmp = this->RunX(input, x, y);
			outData[outX] = ImageUtils::clamp_cast<V>(tmp);
			outX++;
		}

	}

};


template <typename T, typename V>
void ConvolutionFilterSeparable::LoopHelperSingleChannelY(const Image2d<T> & input, Image2d<V> & output,
	int yStart, int yEnd,
	int xStart, int xEnd) const
{
	if (this->kernelDataY.size() < 4)
	{
		//kernel is too small to be vectorized
		ConvolutionFilter::LoopHelper(input, output,
			yStart, yEnd,
			xStart, xEnd,
			*this, &ConvolutionFilterSeparable::RunY);
		return;
	}

	V * outData = output.GetData().data();

#ifdef ENABLE_SIMD

	const T * inData = input.GetData().data();	
	int kn = this->n.GetSize();
	int kn4 = kn - kn % 4;
#endif

	for (int y = yStart; y < yEnd; y++)
	{
		size_t yW = size_t(y) * output.GetWidth();
		size_t outX = size_t(xStart) + yW;

		int xEndSimd = xStart;
#ifdef ENABLE_SIMD

		int len = xEnd - xStart;
		len = len - len % 4;
		xEndSimd = xStart + len;

		size_t y0 = yW - this->n.radius * output.GetWidth();
		size_t y1 = y0 + output.GetWidth();
		size_t y2 = y1 + output.GetWidth();
		size_t y3 = y2 + output.GetWidth();

		for (int x = xStart; x < xEndSimd; x += 4)
		{
			size_t x0 = x + y0;
			size_t x1 = x + y1;
			size_t x2 = x + y2;
			size_t x3 = x + y3;

			__m128 sumDot = _mm_setzero_ps();

			int i = 0;
			for (; i < kn4; i += 4)
			{
				__m128 kx = _mm_set_ps1(this->kernelDataY[i + 0]);
				__m128 ky = _mm_set_ps1(this->kernelDataY[i + 1]);
				__m128 kz = _mm_set_ps1(this->kernelDataY[i + 2]);
				__m128 kw = _mm_set_ps1(this->kernelDataY[i + 3]);


				__m128 dx, dy, dz, dw;

				if constexpr (std::is_same<T, uint8_t>::value)
				{
					//we need co convert uint8_t inputs to float
					__m128i u8_0 = _mm_loadu_si128((const __m128i*)(inData + x0));
					__m128i u8_1 = _mm_loadu_si128((const __m128i*)(inData + x1));
					__m128i u8_2 = _mm_loadu_si128((const __m128i*)(inData + x2));
					__m128i u8_3 = _mm_loadu_si128((const __m128i*)(inData + x3));

					__m128i u32_0 = _mm_unpacklo_epi16(
						_mm_unpacklo_epi8(u8_0, _mm_setzero_si128()),
						_mm_setzero_si128());
					__m128i u32_1 = _mm_unpacklo_epi16(
						_mm_unpacklo_epi8(u8_1, _mm_setzero_si128()),
						_mm_setzero_si128());
					__m128i u32_2 = _mm_unpacklo_epi16(
						_mm_unpacklo_epi8(u8_2, _mm_setzero_si128()),
						_mm_setzero_si128());
					__m128i u32_3 = _mm_unpacklo_epi16(
						_mm_unpacklo_epi8(u8_3, _mm_setzero_si128()),
						_mm_setzero_si128());

					dx = _mm_cvtepi32_ps(u32_0);
					dy = _mm_cvtepi32_ps(u32_1);
					dz = _mm_cvtepi32_ps(u32_2);
					dw = _mm_cvtepi32_ps(u32_3);

				}
				else
				{
					/*
					//load 8 consecutive values
					auto dd = _mm256_loadu_ps(inData + x0);

					//extract parts by shifting and casting to 4 values float
					dx = _mm256_castps256_ps128(dd);
					dy = _mm256_castps256_ps128(_mm256_permutevar8x32_ps(dd, _mm256_set_epi32(0, 0, 0, 0, 4, 3, 2, 1)));
					dz = _mm256_castps256_ps128(_mm256_permutevar8x32_ps(dd, _mm256_set_epi32(0, 0, 0, 0, 5, 4, 3, 2)));
					dw = _mm256_castps256_ps128(_mm256_permutevar8x32_ps(dd, _mm256_set_epi32(0, 0, 0, 0, 6, 5, 4, 3)));
					*/

					dx = _mm_loadu_ps(inData + x0);
					dy = _mm_loadu_ps(inData + x1);
					dz = _mm_loadu_ps(inData + x2);
					dw = _mm_loadu_ps(inData + x3);
				}

				//calculate 4 dots at once
				//[dx, dy, dz, dw] <dot> [kx, ky, kz, kw]

				auto mx = _mm_mul_ps(dx, kx); //dx * kx
				auto my = _mm_fmadd_ps(dy, ky, mx); //mx + dy * ky
				auto mz = _mm_fmadd_ps(dz, kz, my); //my + dz * kz
				auto res = _mm_fmadd_ps(dw, kw, mz); //mz + dw * kw

				sumDot = _mm_add_ps(sumDot, res);

				x0 += 4 * output.GetWidth();
				x1 += 4 * output.GetWidth();
				x2 += 4 * output.GetWidth();
				x3 += 4 * output.GetWidth();
			}

			for (; i < kn; i++)
			{
				auto v = _mm_set_ps1(this->kernelDataX[i]);

				__m128 v2;
				if constexpr (std::is_same<T, uint8_t>::value)
				{
					//we need co convert uint8_t inputs to float
					__m128i u8_0 = _mm_loadu_si128((const __m128i*)(inData + x0));

					__m128i u32_0 = _mm_unpacklo_epi16(
						_mm_unpacklo_epi8(u8_0, _mm_setzero_si128()),
						_mm_setzero_si128());

					v2 = _mm_cvtepi32_ps(u32_0);
				}
				else
				{
					v2 = _mm_loadu_ps(inData + x0);
				}
				sumDot = _mm_add_ps(sumDot, _mm_mul_ps(v, v2));

				x0 += output.GetWidth();
			}

			sumDot = _mm_mul_ps(sumDot, _mm_set_ps1(this->mulX));

			if constexpr (std::is_same<V, uint8_t>::value)
			{
				__m128i asInt = _mm_cvtps_epi32(sumDot);

				asInt = _mm_packus_epi32(asInt, asInt);	// Pack down to 16 bits
				asInt = _mm_packus_epi16(asInt, asInt); // Pack down to 8 bits

				uint32_t res = _mm_cvtsi128_si32(asInt); // Store the lower 32 bits	

				((uint32_t *)(outData + outX))[0] = res;
				outX += 4;
			}
			else
			{
				float tmpRes[4];
				_mm_store_ps(tmpRes, sumDot);

				outData[outX + 0] = tmpRes[0];
				outData[outX + 1] = tmpRes[1];
				outData[outX + 2] = tmpRes[2];
				outData[outX + 3] = tmpRes[3];
				outX += 4;
			}

		}
#endif		

		for (int x = xEndSimd; x < xEnd; x++)
		{
			float tmp = this->RunY(input, x, y);
			outData[outX] = ImageUtils::clamp_cast<V>(tmp);
			outX++;
		}

	}

};

//=================================================================================================


#endif
