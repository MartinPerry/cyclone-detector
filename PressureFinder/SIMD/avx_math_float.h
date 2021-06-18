#ifndef AVX_MATH_FLOAT_H
#define AVX_MATH_FLOAT_H


//https ://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=sqrt&expand=5352,5352&techs=MMX,SSE,SSE2,SSE3,SSSE3,SSE4_1,SSE4_2,AVX,AVX2

#include <immintrin.h>     //AVX2

#include "./avx_utils.h"

	
//=============================================================================
enum class AvxSinCos {
    Sin = 1,
    Cos = 2
};
                                                    
                                                    
template<int Type = int(AvxSinCos::Sin) | int(AvxSinCos::Cos)>
static void _my_mm256_sincos_ps(__m256 x, __m256 *s, __m256 *c)
{
    const float cephes_FOPI = 1.27323954473516f; // 4 / M_PI
    const float minus_cephes_DP1 = -0.78515625f;
    const float minus_cephes_DP2 = -2.4187564849853515625e-4f;
    const float minus_cephes_DP3 = -3.77489497744594108e-8f;

    const float sincof_p0 = -1.9515295891E-4f;
    const float sincof_p1 = 8.3321608736E-3f;
    const float sincof_p2 = -1.6666654611E-1f;
    const float coscof_p0 = 2.443315711809948E-005f;
    const float coscof_p1 = -1.388731625493765E-003f;
    const float coscof_p2 = 4.166664568298827E-002f;
    
	__m256 xmm1, xmm2;
	__m256 sign_bit_sin, y;
    __m256 sign_bit_cos;

	__m256i emm0, emm2;

    __m256i val2 = _mm256_set1_epi32(2);
    __m256i val4 = _mm256_set1_epi32(4);
    __m256i maskAll1 = _mm256_cmpeq_epi32(_mm256_castps_si256(x), _mm256_castps_si256(x)); //fill mask with all 1 => 11111 (32x 1')
    
    if (Type & int(AvxSinCos::Sin))
    {
        // extract the sign bit (upper one)
        //sign_bit_sin = _mm256_and_ps(x, *(__m256*)_i32_256_sign_mask);
        
        // extract the sign bit (upper one)
        __m256i maskSign = _mm256_slli_epi32(maskAll1, 31);     //shift left by 31 bits => 100000... (31x 0')
        sign_bit_sin = _mm256_and_ps(x, _mm256_castsi256_ps(maskSign));
    }
    
	// take the absolute value
	x = _my_mm256_abs_ps(x);
	
	// scale by 4/Pi
	y = _mm256_mul_ps(x, _mm256_set1_ps(cephes_FOPI));

	// store the integer part of y in emm2 (j)
	emm2 = _mm256_cvttps_epi32(y);

	// j=(j+1) & (~1) (see the cephes sources)
    __m256i maskOne = _mm256_srli_epi32(maskAll1, 31);
	emm2 = _mm256_add_epi32(emm2, maskOne);
    emm2 = _mm256_andnot_si256(maskOne, emm2);

	//Convert integers from emm2 in double
	y = _mm256_cvtepi32_ps(emm2);

    
    if (Type & int(AvxSinCos::Cos))
    {
        emm0 = _mm256_sub_epi32(emm2, val2);
        emm0 = _mm256_andnot_si256(emm0, val4);
        emm0 = _mm256_slli_epi32(emm0, 29);
        sign_bit_cos = _mm256_castsi256_ps(emm0);
    }
    
    if (Type & int(AvxSinCos::Sin))
    {
        // get the swap sign flag for the sine
        emm0 = _mm256_and_si256(emm2, val4);
        //left shift by count (29)
        emm0 = _mm256_slli_epi32(emm0, 29);
        __m256 swap_sign_bit_sin = _mm256_castsi256_ps(emm0);
        sign_bit_sin = _mm256_xor_ps(sign_bit_sin, swap_sign_bit_sin);
    }


	// get the polynom selection mask for the sine
	emm2 = _mm256_and_si256(emm2, val2);
	emm2 = _mm256_cmpeq_epi32(emm2, _mm256_setzero_si256());
	__m256 poly_mask = _mm256_castsi256_ps(emm2);


	// The magic pass: "Extended precision modular arithmetic"
	// x = ((x - y * DP1) - y * DP2) - y * DP3;
    x = _mm256_fmadd_ps(y, _mm256_set1_ps(minus_cephes_DP1), x);
    x = _mm256_fmadd_ps(y, _mm256_set1_ps(minus_cephes_DP2), x);
    x = _mm256_fmadd_ps(y, _mm256_set1_ps(minus_cephes_DP3), x);
    

	// Evaluate the first polynom  (0 <= x <= Pi/4)
    __m256 z = _mm256_mul_ps(x, x);
    y = _mm256_set1_ps(coscof_p0);
    y = _mm256_fmadd_ps(y, z, _mm256_set1_ps(coscof_p1));
    y = _mm256_fmadd_ps(y, z, _mm256_set1_ps(coscof_p2));
    
	y = _mm256_mul_ps(y, z);
	y = _mm256_mul_ps(y, z);
    
    //create constant 0.5 (see Agner Fog)
    __m256i tmpi = _mm256_slli_epi32(maskAll1, 26);
	__m256 val0p5 = _mm256_castsi256_ps(_mm256_srli_epi32(tmpi, 2));
    
    //create constant 1.0 (see Agner Fog)
    tmpi = _mm256_slli_epi32(maskAll1, 25);
	__m256 val1p0 = _mm256_castsi256_ps(_mm256_srli_epi32(tmpi, 2));
    
    //fnmadd = -(a*b)+c
    y = _mm256_fnmadd_ps(z, val0p5, y);
	//__m256 tmp = _mm256_mul_ps(z, val0p5);
	//y = _mm256_sub_ps(y, tmp);
	y = _mm256_add_ps(y, val1p0);

	// Evaluate the second polynom  (Pi/4 <= x <= 0)
    __m256 y2 = _mm256_set1_ps(sincof_p0);
    y2 = _mm256_fmadd_ps(y2, z, _mm256_set1_ps(sincof_p1));
    y2 = _mm256_fmadd_ps(y2, z, _mm256_set1_ps(sincof_p2));
    y2 = _mm256_mul_ps(y2, z);
    y2 = _mm256_fmadd_ps(y2, x, x);

    
	// select the correct result from the two polynoms
	__m256 ysin2 = _mm256_and_ps(poly_mask, y2);
	__m256 ysin1 = _mm256_andnot_ps(poly_mask, y);

	// update the sign
    if (Type & int(AvxSinCos::Sin))
    {
        xmm1 = _mm256_add_ps(ysin1, ysin2);
        *s = _mm256_xor_ps(xmm1, sign_bit_sin);
    }
    if (Type & int(AvxSinCos::Cos))
    {
        y2 = _mm256_sub_ps(y2, ysin2);
        y = _mm256_sub_ps(y, ysin1);
        xmm2 = _mm256_add_ps(y, y2);
        *c = _mm256_xor_ps(xmm2, sign_bit_cos);
    }
}

static __m256 _my_mm256_sin_ps(__m256 x)
{
	__m256 s, c;
    _my_mm256_sincos_ps<int(AvxSinCos::Sin)>(x, &s, &c);
	return s;
}

static __m256 _my_mm256_cos_ps(__m256 x)
{
	__m256 s, c;
	_my_mm256_sincos_ps<int(AvxSinCos::Cos)>(x, &s, &c);
	return c;
}
                                                     
static __m256 _my_mm256_tan_ps(__m256 x)
{
    __m256 s, c;
    _my_mm256_sincos_ps(x, &s, &c);
    return _mm256_div_ps(s, c);
}

//=============================================================================

template<int Type = int(AvxSinCos::Sin) | int(AvxSinCos::Cos)>
static void _my_mm256_asincos_ps(__m256 x,  __m256 *s, __m256 *c)
{
    const __m256 PIO2F = _mm256_set1_ps(1.5707963267948966192f);
    const __m256 c_05 = _mm256_set1_ps(0.5);
    const __m256 c_1 = _mm256_set1_ps(1.0);
    
    //__m256i t1 = _mm256_castps_si256(x);    // reinterpret as 32-bit integer
    //__m256i t2 = _mm256_srli_epi32(t1, 31);            // extend sign bit
    //__m256 signBit = _mm256_castsi256_ps(t2);       // reinterpret as 32-bit Boolean
    
    
    __m256 signBit = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LT_OS);
        
    x = _my_mm256_abs_ps(x);
    
    //test if a is > 0.5
    __m256 over05 = _mm256_cmp_ps(x, c_05, _CMP_GT_OS);
        
    __m256 z1 = _mm256_sub_ps(c_1, x);
    z1 = _mm256_mul_ps(z1, c_05);
    __m256 x1 = _mm256_sqrt_ps(z1);
        
    __m256 x2 = x;
    __m256 z2 = _mm256_mul_ps(x2, x2);
        
    x = _my_mm256_select(over05, x1, x2);
    __m256 z = _my_mm256_select(over05, z1, z2);
    
    __m256 pz = _mm256_fmadd_ps(z, _mm256_set1_ps(4.2163199048E-2f), _mm256_set1_ps(2.4181311049E-2f));
    pz = _mm256_fmadd_ps(pz, z, _mm256_set1_ps(4.5470025998E-2f));
    pz = _mm256_fmadd_ps(pz, z, _mm256_set1_ps(7.4953002686E-2f));
    pz = _mm256_fmadd_ps(pz, z, _mm256_set1_ps(1.6666752422E-1f));
    pz = _mm256_mul_ps(pz, z);
    z = _mm256_fmadd_ps(pz, x, x);
    
    //__m256 pz = _mm256_mul_ps(z, _mm256_set1_ps(4.2163199048E-2));
    //pz = _mm256_add_ps(pz, _mm256_set1_ps(2.4181311049E-2));
    //pz = _mm256_mul_ps(pz, z);
    //pz = _mm256_add_ps(pz, _mm256_set1_ps(4.5470025998E-2));
    //pz = _mm256_mul_ps(pz, z);
    //pz = _mm256_add_ps(pz, _mm256_set1_ps(7.4953002686E-2));
    //pz = _mm256_mul_ps(pz, z);
    //pz = _mm256_add_ps(pz, _mm256_set1_ps(1.6666752422E-1));
    //pz = _mm256_mul_ps(pz, z);
    //pz = _mm256_mul_ps(pz, x);
    //z = _mm256_add_ps(pz, x);
        
    __m256 tmp2z = _mm256_add_ps(z, z);
    
    if (Type & int(AvxSinCos::Cos))
    {
        const __m256 PIF = _mm256_set1_ps(3.14159265358979323846f);
        
        __m256 tmp = _my_mm256_select(signBit, _mm256_sub_ps(PIF, tmp2z), tmp2z);
        __m256 tmp2 = _mm256_sub_ps(PIO2F, _my_mm256_select(signBit, _my_mm256_swap_sign(z), z));
        *c = _my_mm256_select(over05, tmp, tmp2);
    }
    
    if (Type & int(AvxSinCos::Sin))
    {
        __m256 tmp = _mm256_sub_ps(PIO2F, tmp2z);
        z1 = _my_mm256_select(over05, tmp, z);
        *s = _my_mm256_select(signBit, _my_mm256_swap_sign(z1), z1);
    }
}

static __m256 _my_mm256_asin_ps(__m256 x)
{
    __m256 s, c;
    _my_mm256_asincos_ps<int(AvxSinCos::Sin)>(x, &s, &c);
    return s;
}
                                       
static __m256 _my_mm256_acos_ps(__m256 x)
{
    __m256 s, c;
    _my_mm256_asincos_ps<int(AvxSinCos::Cos)>(x, &s, &c);
    return c;
}
                                       
static __m256 _my_mm256_atan_ps(__m256 x)
{
    const float PIO2F = 1.5707963267948966192f;
    const float PIO4F = 0.7853981633974483096f;
    const float TAN_3_PI_DIV_8 = 2.414213562373095f;
    const float TAN_PI_DIV_8 = 0.4142135623730950f;
    
    const __m256 c_m1 = _mm256_set1_ps(-1.0);
    
    //__m256i t1 = _mm256_castps_si256(x);    // reinterpret as 32-bit integer
    //__m256i t2 = _mm256_srli_epi32(t1, 31);            // extend sign bit
    //__m256 signBit = _mm256_castsi256_ps(t2);       // reinterpret as 32-bit Boolean
    
    
    __m256 signBit = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LT_OS);
    
    x = _my_mm256_abs_ps(x);
    
    // small:  x < TAN_PI_DIV_8
    // medium: TAN_PI_DIV_8 <= x <= TAN_3_PI_DIV_8
    // big:    x > TAN_3_PI_DIV_8
    __m256 big = _mm256_cmp_ps(x, _mm256_set1_ps(TAN_3_PI_DIV_8), _CMP_GT_OS);
    __m256 medium = _mm256_cmp_ps(x, _mm256_set1_ps(TAN_PI_DIV_8), _CMP_GT_OS);
    
    __m256 s = _my_mm256_select(big, _mm256_set1_ps(PIO2F),
		_my_mm256_select(medium, _mm256_set1_ps(PIO4F), _mm256_setzero_ps()));
    
    __m256 tDiv = _mm256_div_ps(c_m1, x); //-1./x
    
    __m256 t1 = _mm256_add_ps(x, c_m1); //x+(-1) => x-1
    __m256 t2 = _mm256_sub_ps(x, c_m1); //x-(-1) => x+1
    
    x = _my_mm256_select(big, tDiv, _my_mm256_select(medium, _mm256_div_ps(t1, t2), x));
    
    __m256 z =_mm256_mul_ps(x, x);
    
    //s += ((( 8.05374449538e-2 * z - 1.38776856032E-1) * z + 1.99777106478E-1) * z - 3.33329491539E-1) * z * x + x;
    
    __m256 pz = _mm256_fmadd_ps(z, _mm256_set1_ps(8.05374449538e-2f), _mm256_set1_ps(-1.38776856032E-1f));
    pz = _mm256_fmadd_ps(pz, z, _mm256_set1_ps(1.99777106478E-1f));
    pz = _mm256_fmadd_ps(pz, z, _mm256_set1_ps(-3.33329491539E-1f));
    pz = _mm256_mul_ps(pz, z);
    z = _mm256_fmadd_ps(pz, x, x);
    s = _mm256_add_ps(s, z);
    
    return _my_mm256_select(signBit, _my_mm256_swap_sign(s), s);
}

static __m256 _my_mm256_atan2_ps(__m256 y, __m256 x)
{
	const float PIF = 3.14159265358979323846f;
	const float PIO2F = 1.5707963267948966192f;
	const float PIO4F = 0.7853981633974483096f;	
	const float TAN_PI_DIV_8 = 0.4142135623730950f;

	const __m256 c_m1 = _mm256_set1_ps(-1.0);

	__m256 signBitX = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LT_OS);
	__m256 signBitY = _mm256_cmp_ps(y, _mm256_setzero_ps(), _CMP_LT_OS);

	__m256 x1 = _my_mm256_abs_ps(x);
	__m256 y1 = _my_mm256_abs_ps(y);

	__m256 swapxy = _mm256_cmp_ps(y1, x1, _CMP_GT_OS);
	__m256 x2 = _my_mm256_select(swapxy, y1, x1);
	__m256 y2 = _my_mm256_select(swapxy, x1, y1);
	__m256 t = _mm256_div_ps(y2, x2);

	__m256 medium = _mm256_cmp_ps(t, _mm256_set1_ps(TAN_PI_DIV_8), _CMP_GT_OS);
	__m256 s = _my_mm256_select(medium, _mm256_set1_ps(PIO4F), _mm256_setzero_ps());
	__m256 t1 = _mm256_add_ps(t, c_m1); //t+(-1) => t-1
	__m256 t2 = _mm256_sub_ps(t, c_m1); //t-(-1) => t+1
	__m256 z = _my_mm256_select(medium, _mm256_div_ps(t1, t2), t);
	__m256 zz = _mm256_mul_ps(z, z);

	__m256 pz = _mm256_fmadd_ps(zz, _mm256_set1_ps(8.05374449538e-2f), _mm256_set1_ps(-1.38776856032E-1f));
	pz = _mm256_fmadd_ps(pz, zz, _mm256_set1_ps(1.99777106478E-1f));
	pz = _mm256_fmadd_ps(pz, zz, _mm256_set1_ps(-3.33329491539E-1f));
	pz = _mm256_mul_ps(pz, zz);
	z = _mm256_fmadd_ps(pz, z, z);
	s = _mm256_add_ps(s, z);

	s = _my_mm256_select(swapxy, _mm256_sub_ps(_mm256_set1_ps(PIO2F), s), s);
	s = _my_mm256_select(signBitX, _mm256_sub_ps(_mm256_set1_ps(PIF), s), s);

	//re = select((x | y) == 0.f, 0.f, re);    // atan2(0,0) = 0 by convention

	return _my_mm256_select(signBitY, _my_mm256_swap_sign(s), s);
}
                                       
//=============================================================================
                                       
static __m256 _my_mm256_exp_ps(__m256 x)
{
    const float exp_hi = 88.3762626647949f;
    const float exp_lo = -88.3762626647949f;
    const float cephes_LOG2EF = 1.44269504088896341f;
    const float cephes_exp_C1 = 0.693359375f;
    const float cephes_exp_C2 = -2.12194440e-4f;
    const float cephes_exp_p0 = 1.9875691500E-4f;
    const float cephes_exp_p1 = 1.3981999507E-3f;
    const float cephes_exp_p2 = 8.3334519073E-3f;
    const float cephes_exp_p3 = 4.1665795894E-2f;
    const float cephes_exp_p4 = 1.6666665459E-1f;
    const float cephes_exp_p5 = 5.0000001201E-1f;
    
	__m256i emm0;

    __m256i maskAll1 = _mm256_cmpeq_epi32(_mm256_castps_si256(x), _mm256_castps_si256(x)); //fill mask with all 1 => 11111 (32x 1')
    
    //create constant 0.5 (see agner Fog)
    __m256i tmpi = _mm256_slli_epi32(maskAll1, 26);
	__m256 val0p5 = _mm256_castsi256_ps(_mm256_srli_epi32(tmpi, 2));
    
    //create constant 1.0 (see Agner Fog)
	tmpi = _mm256_slli_epi32(maskAll1, 25);
    __m256 val1p0 = _mm256_castsi256_ps(_mm256_srli_epi32(tmpi, 2));
    
    //=====================
    
	x = _mm256_min_ps(x, _mm256_set1_ps(exp_hi));
	x = _mm256_max_ps(x, _mm256_set1_ps(exp_lo));

	// express exp(x) as exp(g + n*log(2))
    __m256 fx = _mm256_fmadd_ps(x, _mm256_set1_ps(cephes_LOG2EF), val0p5);

	fx = _mm256_floor_ps(fx);

    //=======
	// how to perform a floorf with SSE: just below
	//emm0 = _mm256_cvttpd_epi32(fx);
	//tmp = _mm256_cvtepi32_pd(emm0);

	// if greater, substract 1
	//__m256d mask = _mm_cmpgt_pd(tmp, fx);
	//mask = _mm256_and_pd(mask, one);
	//fx = _mm256_sub_pd(tmp, mask);
    //=======
    
    //fnmadd = -(a*b)+c
    x = _mm256_fnmadd_ps(fx, _mm256_set1_ps(cephes_exp_C1), x);
    x = _mm256_fnmadd_ps(fx, _mm256_set1_ps(cephes_exp_C2), x);

	__m256 z = _mm256_mul_ps(x, x);

    __m256 y = _mm256_fmadd_ps(_mm256_set1_ps(cephes_exp_p0), x, _mm256_set1_ps(cephes_exp_p1));
    y = _mm256_fmadd_ps(y, x, _mm256_set1_ps(cephes_exp_p2));
    y = _mm256_fmadd_ps(y, x, _mm256_set1_ps(cephes_exp_p3));
    y = _mm256_fmadd_ps(y, x, _mm256_set1_ps(cephes_exp_p4));
    y = _mm256_fmadd_ps(y, x, _mm256_set1_ps(cephes_exp_p5));
    y = _mm256_fmadd_ps(y, z, x);
    y = _mm256_add_ps(y, val1p0);
    
    /*
	__m256 y = _mm256_set1_ps(cephes_exp_p0);
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_exp_p1));
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_exp_p2));
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_exp_p3));
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_exp_p4));
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_exp_p5));
	y = _mm256_mul_ps(y, z);
	y = _mm256_add_ps(y, x);
	y = _mm256_add_ps(y, val1p0);
     */
    
	// build 2^n
	emm0 = _mm256_cvttps_epi32(fx);
	emm0 = _mm256_add_epi32(emm0, _mm256_set1_epi32(0x7f));
	emm0 = _mm256_slli_epi32(emm0, 23);
	__m256 pow2n = _mm256_castsi256_ps(emm0);

	y = _mm256_mul_ps(y, pow2n);

	return y;

}


static __m256 _my_mm256_log_ps(__m256 x) {

	const int MANTISA_SIZE = 23;

    const float cephes_SQRTHF = 0.707106781186547524f;
    const float cephes_log_p0 = 7.0376836292E-2f;
    const float cephes_log_p1 = -1.1514610310E-1f;
    const float cephes_log_p2 = 1.1676998740E-1f;
    const float cephes_log_p3 = -1.2420140846E-1f;
    const float cephes_log_p4 = +1.4249322787E-1f;
    const float cephes_log_p5 = -1.6668057665E-1f;
    const float cephes_log_p6 = +2.0000714765E-1f;
    const float cephes_log_p7 = -2.4999993993E-1f;
    const float cephes_log_p8 = +3.3333331174E-1f;
    const float cephes_log_q1 = -2.12194440e-4f;
    const float cephes_log_q2 = 0.693359375f;
    
	__m256i emm0;
    
    __m256i maskAll1 = _mm256_cmpeq_epi32(_mm256_castps_si256(x), _mm256_castps_si256(x)); //fill mask with all 1 => 11111 (32x 1')
    
    //create constant 0.5 (see agner Fog)
    __m256i tmpi = _mm256_slli_epi32(maskAll1, 26);
	__m256 val0p5 = _mm256_castsi256_ps(_mm256_srli_epi32(tmpi, 2));
    
    //create constant 1.0 (see Agner Fog)
	tmpi = _mm256_slli_epi32(maskAll1, 25);
	__m256 val1p0 = _mm256_castsi256_ps(_mm256_srli_epi32(tmpi, 2));
    
    //=====================

	__m256 invalid_mask = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LE_OQ);
	
	x = _mm256_max_ps(x, _mm256_castsi256_ps(_mm256_set1_epi32(0x7f)));  /* cut off denormalized stuff */

	//__m256i tmpi = _mm256_castps_si256(x);

	emm0 = _mm256_srli_epi32(_mm256_castps_si256(x), MANTISA_SIZE);

	/* keep only the fractional part */
	x = _mm256_and_ps(x, _mm256_castsi256_ps(_mm256_set1_epi32(static_cast<uint32_t>(~0x7f80'0000))));
	x = _mm256_or_ps(x,  val0p5);


	emm0 = _mm256_sub_epi32(emm0, _mm256_set1_epi32(0x7f));
	__m256 e = _mm256_cvtepi32_ps(emm0);

	e = _mm256_add_ps(e, val1p0);

	/* part2:
	   if( x < SQRTHF ) {
		 e -= 1;
		 x = x + x - 1.0;
	   } else { x = x - 1.0; }
	*/
	__m256 mask = _mm256_cmp_ps(x, _mm256_set1_ps(cephes_SQRTHF), _CMP_LT_OQ);
	__m256 tmp = _mm256_and_ps(x, mask);
	x = _mm256_sub_ps(x, val1p0);
	e = _mm256_sub_ps(e, _mm256_and_ps(val1p0, mask));
	x = _mm256_add_ps(x, tmp);


	__m256 z = _mm256_mul_ps(x, x);

    __m256 y = _mm256_fmadd_ps(_mm256_set1_ps(cephes_log_p0), x, _mm256_set1_ps(cephes_log_p1));
    y = _mm256_fmadd_ps(y, x, _mm256_set1_ps(cephes_log_p2));
    y = _mm256_fmadd_ps(y, x, _mm256_set1_ps(cephes_log_p3));
    y = _mm256_fmadd_ps(y, x, _mm256_set1_ps(cephes_log_p4));
    y = _mm256_fmadd_ps(y, x, _mm256_set1_ps(cephes_log_p5));
    y = _mm256_fmadd_ps(y, x, _mm256_set1_ps(cephes_log_p6));
    y = _mm256_fmadd_ps(y, x, _mm256_set1_ps(cephes_log_p7));
    y = _mm256_fmadd_ps(y, x, _mm256_set1_ps(cephes_log_p8));
    y = _mm256_mul_ps(y, x);
    y = _mm256_mul_ps(y, z);
    
    y = _mm256_fmadd_ps(e, _mm256_set1_ps(cephes_log_q1), y);
    y = _mm256_fnmadd_ps(z, val0p5, y);
    
    
    /*
	__m256 y = _mm256_set1_ps(cephes_log_p0);
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_log_p1));
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_log_p2));
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_log_p3));
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_log_p4));
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_log_p5));
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_log_p6));
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_log_p7));
	y = _mm256_mul_ps(y, x);
	y = _mm256_add_ps(y, _mm256_set1_ps(cephes_log_p8));
	y = _mm256_mul_ps(y, x);

	y = _mm256_mul_ps(y, z);
    

	tmp = _mm256_mul_ps(e, _mm256_set1_ps(cephes_log_q1));
	y = _mm256_add_ps(y, tmp);


	tmp = _mm256_mul_ps(z, val0p5);
	y = _mm256_sub_ps(y, tmp);
    */

	x = _mm256_add_ps(x, y);
    x = _mm256_fmadd_ps(e, _mm256_set1_ps(cephes_log_q2), x);
    //tmp = _mm256_mul_ps(e, _mm256_set1_ps(cephes_log_q2));
	//x = _mm256_add_ps(x, tmp);
	x = _mm256_or_ps(x, invalid_mask); // negative arg will be NAN
	return x;
}

static __m256 _my_mm256_pow_ps(const __m256 & x, const __m256 & y)
{
    __m256 tmp = _my_mm256_log_ps(x);
    return _my_mm256_exp_ps(_mm256_mul_ps(y, tmp));
}

#endif