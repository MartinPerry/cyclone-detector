#ifndef AVX_UTILS_H
#define AVX_UTILS_H

#include <cstdint>
#include <algorithm>

#include <immintrin.h>     //AVX2

static const int MM256_ELEMENT_COUNT = 8;
static const int MM256d_ELEMENT_COUNT = 4;

//it is with define in immintrin.h on windows
//was not compiling on Linux without this (why?)
#ifndef _mm256_cvtsi256_si32
#	define _mm256_cvtsi256_si32(a) (_mm_cvtsi128_si32(_mm256_castsi256_si128(a)))
#endif

//=====================================================================================
// double based operations
//=====================================================================================


/// <summary>
/// Convert input 4-double values to 4-uint8_t values with 
/// clamping of values outside range to [0, 255]
/// They are packed inside uin32_t
/// For conversion, rounding mode is used. 
/// Rounding mode is set via _MM_SET_ROUNDING_MODE
/// Example for round: _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
/// 
/// </summary>
/// <param name="v"></param>
/// <returns></returns>
static inline uint32_t _my_mm256d_cast_clamp_to_uint8(const __m256d & v)
{	
	__m128i asInt = _mm256_cvtpd_epi32(v);
	
	asInt = _mm_packus_epi32(asInt, asInt);	// Pack down to 16 bits
	asInt = _mm_packus_epi16(asInt, asInt); // Pack down to 8 bits

	return _mm_cvtsi128_si32(asInt); // Store the lower 32 bits		
};

//=====================================================================================
// float based operations
//=====================================================================================

///<summary>
/// Select and return a or b based on mask value.
/// If mask value is 1, return a; else return b
///</summary>
static inline __m256 _my_mm256_select(const __m256 & mask, const __m256 & a, const __m256 & b)
{
	return  _mm256_blendv_ps(b, a, mask);
}

///<summary>
/// Select and return a or b based on mask value.
/// If v1 COMPARISON_OPERATOR v2 return a; else return b
/// COMPARISON_OPERATOR - https://stackoverflow.com/questions/16988199/how-to-choose-avx-compare-predicate-variants
///</summary>
template <int COMPARISON_OPERATOR>
static inline __m256 _my_mm256_select_if(const __m256 & v1, const __m256 & v2,
	const __m256 & a, const __m256 & b)
{
	__m256 mask = _mm256_cmp_ps(v1, v2, COMPARISON_OPERATOR);
	return _my_mm256_select(mask, a, b);
};


/// <summary>
/// Convert input 8-float values to 8-uint8_t values with 
/// clamping of values outside range to [0, 255]
/// They are packed inside uin64_t
/// For conversion, rounding mode is used. 
/// Rounding mode is set via _MM_SET_ROUNDING_MODE
/// Example for round: _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
/// 
/// </summary>
/// <param name="v"></param>
/// <returns></returns>
static inline uint64_t _my_mm256_cast_clamp_to_uint8(const __m256 & v)
{
	__m256i asInt = _mm256_cvtps_epi32(v);

	asInt = _mm256_packus_epi32(asInt, asInt); // Pack down to 16 bits
	asInt = _mm256_packus_epi16(asInt, asInt); // Pack down to 8 bits

	uint32_t v0 = _mm256_cvtsi256_si32(asInt); // Store the lower 32 bits		

	asInt = _mm256_permute2f128_si256(asInt, asInt, 1); //swap lower and high 32bits

	uint32_t v1 = _mm256_cvtsi256_si32(asInt); // Store the upper 32 bits		

	return (uint64_t(v1) << 32) | uint64_t(v0);
};



/// <summary>
/// Calculate abs of 8 float
/// and get signBit of input 
/// It can be used to "unabs" result back via:
/// _mm256_xor_ps(absV, signBit);
/// -> we put negative sign back
/// Example:
/// v = (-5, 5, -4) => res = (5, 5, 4)
/// xor(res, signBit) = (-5, 5, -4)
/// </summary>
/// <param name="v"></param>
/// <param name="signBit"></param>
/// <returns></returns>
static inline __m256 _my_mm256_abs_ps(const __m256 & v, __m256 & signBit)
{
	__m256i maskAll1 = _mm256_cmpeq_epi32(_mm256_castps_si256(v), _mm256_castps_si256(v)); //fill mask with all 1 => 11111 (32x 1')	
	__m256i maskSign = _mm256_slli_epi32(maskAll1, 31);  //shift left by 31 bits => 100000... (31x 0')
	signBit = _mm256_and_ps(v, _mm256_castsi256_ps(maskSign));
	
	return _mm256_xor_ps(v, signBit); //remove negative sign bit
}


/// <summary>
/// Calculate abs of 8 float
/// with bitmasking sign bit
/// </summary>
static inline __m256 _my_mm256_abs_ps(const __m256 & v)
{
	__m256i mask = _mm256_cmpeq_epi32(_mm256_castps_si256(v), _mm256_castps_si256(v)); //fill mask with all 1 => 11111 (32x 1')
	mask = _mm256_srli_epi32(mask, 1);     //shift right by 1 bits => 011111... (31x 1')
	return _mm256_and_ps(v, _mm256_castsi256_ps(mask)); //and to "remove" sign bit
	//return _mm256_and_ps(v, *(__m256*)_i32_256_inv_sign_mask); //0x7FFF'FFFF
}

static inline __m256 _my_mm256_swap_sign(const __m256 & v)
{
	__m256i mask = _mm256_cmpeq_epi32(_mm256_castps_si256(v), _mm256_castps_si256(v)); //fill mask with all 1 => 11111 (32x 1')
	mask = _mm256_slli_epi32(mask, 31);     //shift left by 31 bits => 100000... (31x 0')
	return _mm256_xor_ps(v, _mm256_castsi256_ps(mask)); //xor to change sign bit

	//uint32_t c = static_cast<uint32_t>(~0x7FFF'FFFF);
	//static const ALIGN(32) uint32_t _i32_sign_mask[8] = { c, c, c, c, c, c, c, c };
	//return _mm256_xor_ps(v, *(__m256*)_i32_sign_mask); //~0x7FFF'FFFF
}


/// <summary>
/// clamp v to interval [min, max]
/// uses min() and max() functions to do this
/// </summary>
/// <param name="v"></param>
/// <param name="min"></param>
/// <param name="max"></param>
/// <returns></returns>
static inline __m256 _my_mm256_clamp(const __m256 & v, const __m256 & min, const __m256 & max)
{
	return _mm256_max_ps(min, _mm256_min_ps(v, max));;
}

/// <summary>
/// Simple helper method to get horizontal max
/// it will unpack vector to array and iterate array
/// </summary>
/// <param name="v"></param>
/// <returns></returns>
static inline float _my_mm256_simple_hmax(const __m256 & v)
{
	float vals[MM256_ELEMENT_COUNT];
	_mm256_storeu_ps(vals, v);

	float r = vals[0];
	for (int i = 1; i < MM256_ELEMENT_COUNT; i++)
	{
		r = std::max(r, vals[i]);
	}

	return r;
}

/// <summary>
/// Simple helper method to get horizontal min
/// it will unpack vector to array and iterate array
/// </summary>
/// <param name="v"></param>
/// <returns></returns>
static inline float _my_mm256_simple_hmin(const __m256 & v)
{
	float vals[MM256_ELEMENT_COUNT];
	_mm256_storeu_ps(vals, v);

	float r = vals[0];
	for (int i = 1; i < MM256_ELEMENT_COUNT; i++)
	{
		r = std::min(r, vals[i]);
	}

	return r;
}

//=====================================================================================
// uint8_t based operations
// mostly based on http://www.alfredklomp.com/programming/sse-intrinsics/
//=====================================================================================

/// <summary>
/// Returns 0xFF where x > y:
/// http://www.alfredklomp.com/programming/sse-intrinsics/
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
static inline __m128i _my_mm_cmpgt_epu8(__m128i x, __m128i y)
{
	// Returns 0xFF where x > y:
	return _mm_andnot_si128(
		_mm_cmpeq_epi8(x, y),
		_mm_cmpeq_epi8(_mm_max_epu8(x, y), x)
	);
}

/// <summary>
/// Returns 0xFF where x < y:
/// http://www.alfredklomp.com/programming/sse-intrinsics/
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
static inline __m128i _my_mm_cmplt_epu8(__m128i x, __m128i y)
{
	// Returns 0xFF where x < y:
	return _my_mm_cmpgt_epu8(y, x);
}

/// <summary>
/// Returns 0xFF where x <= y:
/// http://www.alfredklomp.com/programming/sse-intrinsics/
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
static inline __m128i _my_mm_cmple_epu8(__m128i x, __m128i y)
{
	// Returns 0xFF where x <= y:
	return _mm_cmpeq_epi8(_mm_min_epu8(x, y), x);
}

/// <summary>
/// Returns 0xFF where x >= y:
/// http://www.alfredklomp.com/programming/sse-intrinsics/
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
static inline __m128i _my_mm_cmpge_epu8(__m128i x, __m128i y)
{
	// Returns 0xFF where x >= y:
	return _my_mm_cmple_epu8(y, x);
}



/*
//https://stackoverflow.com/questions/8193601/sse-multiplication-16-x-uint8-t
__m128i mullo_epi8(__m128i a, __m128i b)
{
	// unpack and multiply
	__m128i dst_even = _mm_mullo_epi16(a, b);
	__m128i dst_odd = _mm_mullo_epi16(_mm_srli_epi16(a, 8), _mm_srli_epi16(b, 8));
	// repack
#ifdef __AVX2__
	// only faster if have access to VPBROADCASTW
		return _mm_or_si128(_mm_slli_epi16(dst_odd, 8), _mm_and_si128(dst_even, _mm_set1_epi16(0xFF)));
#else
		return _mm_or_si128(_mm_slli_epi16(dst_odd, 8), _mm_srli_epi16(_mm_slli_epi16(dst_even, 8), 8));
#endif
	}
	*/

#endif
