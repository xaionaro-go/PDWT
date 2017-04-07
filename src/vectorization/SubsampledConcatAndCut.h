#ifndef SUBSAMPLEDCONCATANDCUT_H
#define SUBSAMPLEDCONCATANDCUT_H

// Local
#include "MetaHelper.h"
#include "vectorization.h"

template<typename T, class VecT, int SHIFT>
struct SubsampledConcatAndCut {
  static VecT  Concat( VecT a, VecT b, VecT c) {
    assert(false);
  }
  static VecT  Concat( VecT a, VecT b) {
    assert(false);
  }
};

#ifdef USE_AVX
template<>
struct SubsampledConcatAndCut<float,__m128,0> {
  static __m128  Concat( __m128 a, __m128 b, __m128 c) {
    return _mm_shuffle_ps(a,b,0b10001000);
  }
};
template<>
struct SubsampledConcatAndCut<float,__m128,1> {
  static __m128  Concat( __m128 a, __m128 b, __m128 c) {
    return _mm_shuffle_ps(a,b,0b11011101);
  }
};
template<>
struct SubsampledConcatAndCut<float,__m128,2> {
  static __m128  Concat( __m128 a, __m128 b, __m128 c) {
    auto x = _mm_permute_ps(_mm_shuffle_ps(a,b,0b10000010),0b01111000);
    return _mm_blend_ps(x, _mm_permute_ps(c,57),0b00001000);
  }
};
template<>
struct SubsampledConcatAndCut<float,__m128,3> {
  static __m128  Concat( __m128 a, __m128 b, __m128 c) {
    auto x = _mm_permute_ps(_mm_shuffle_ps(a,b,0b11010011),0b01111000);
    return _mm_blend_ps(x, _mm_permute_ps(c,120),0b00001000);
  }
};
template<>
struct SubsampledConcatAndCut<double,__m128d,0> {
  static __m128d  Concat( __m128d a, __m128d b) {
    return _mm_shuffle_pd(a,b,0b00000000);
  }
};
template<>
struct SubsampledConcatAndCut<double,__m128d,1> {
  static __m128d  Concat( __m128d a, __m128d b) {
    return _mm_shuffle_pd(a,b,0b00001111);
  }
};
#elif defined USE_AVX2
template<>
struct SubsampledConcatAndCut<float,__m256,0> {
  static __m256  Concat( __m256 a, __m256 b, __m256 c) {
    return (__m256)_mm256_permute4x64_epi64((__m256i)
      _mm256_shuffle_ps(a,b,0b10001000),216);
  }
};
template<>
struct SubsampledConcatAndCut<float,__m256,1> {
  static __m256  Concat( __m256 a, __m256 b, __m256 c) {
    return (__m256)_mm256_permute4x64_epi64((__m256i)
      _mm256_shuffle_ps(a,b,0b11011101),216);
  }
};
template<>
struct SubsampledConcatAndCut<float,__m256,2> {
  static __m256  Concat( __m256 a, __m256 b, __m256 c) {
    a=_mm256_permutevar8x32_ps(a,
      _mm256_set_epi32(0,0,0,0,0,6,4,2));
    b=_mm256_permutevar8x32_ps(b,
      _mm256_set_epi32(0,6,4,2,0,0,0,0));
    a=_mm256_blend_ps(a,b,0b01111000);
    return _mm256_blend_ps(a,_mm256_permutevar8x32_ps(c,
      _mm256_set_epi32(0,0,0,0,0,0,0,0)),0b10000000);
  }
};
template<>
struct SubsampledConcatAndCut<float,__m256,3> {
  static __m256  Concat( __m256 a, __m256 b, __m256 c) {
    a=_mm256_permutevar8x32_ps(a,
      _mm256_set_epi32(0,0,0,0,0,7,5,3));
    b=_mm256_permutevar8x32_ps(b,
      _mm256_set_epi32(0,7,5,3,1,0,0,0));
    a=_mm256_blend_ps(a,b,0b01111000);
    return _mm256_blend_ps(a,_mm256_permutevar8x32_ps(c,
      _mm256_set_epi32(1,0,0,0,0,0,0,0)),0b10000000);
  }
};
template<>
struct SubsampledConcatAndCut<float,__m256,4> {
  static __m256  Concat( __m256 a, __m256 b, __m256 c) {
    auto x = (__m256) _mm256_permute2x128_si256(
      (__m256i)_mm256_permute_ps(a,0b11011000),
      (__m256i)_mm256_permute_ps(c,0b10001101),97);
    auto y = (__m256)_mm256_permute4x64_epi64((__m256i)
      _mm256_permute_ps(b,0b10001101),180);
    return _mm256_blend_ps(x,y,0b00111100);
  }
};
template<>
struct SubsampledConcatAndCut<float,__m256,5> {
  static __m256  Concat( __m256 a, __m256 b, __m256 c) {
    auto x = (__m256) _mm256_permute2x128_si256(
      (__m256i)_mm256_permute_ps(a,0b10001101),
      (__m256i)_mm256_permute_ps(c,0b11011000),97);
    auto y = (__m256)_mm256_permute4x64_epi64((__m256i)
      _mm256_permute_ps(b,0b11011000),180);
    return _mm256_blend_ps(x,y,0b00111100);
  }
};
template<>
struct SubsampledConcatAndCut<float,__m256,6> {
  static __m256  Concat( __m256 a, __m256 b, __m256 c) {
    a=_mm256_permutevar8x32_ps(a,
      _mm256_set_epi32(0,0,0,0,0,0,0,6));
    b=_mm256_permutevar8x32_ps(b,
      _mm256_set_epi32(0,0,0,6,4,2,0,0));
    a=_mm256_blend_ps(a,b,0b00011110);
    return _mm256_blend_ps(a,_mm256_permutevar8x32_ps(c,
      _mm256_set_epi32(4,2,0,0,0,0,0,0)),0b11100000);
  }
};
template<>
struct SubsampledConcatAndCut<float,__m256,7> {
  static __m256  Concat( __m256 a, __m256 b, __m256 c) {
    a=_mm256_permutevar8x32_ps(a,
      _mm256_set_epi32(0,0,0,0,0,0,0,7));
    b=_mm256_permutevar8x32_ps(b,
      _mm256_set_epi32(0,0,0,7,5,3,1,0));
    a=_mm256_blend_ps(a,b,0b00011110);
    return _mm256_blend_ps(a,_mm256_permutevar8x32_ps(c,
      _mm256_set_epi32(5,3,1,0,0,0,0,0)),0b11100000);
  }
};
template<>
struct SubsampledConcatAndCut<double,__m256d,0> {
  static __m256d  Concat( __m256d a, __m256d b, __m256d c) {
    return (__m256d) _mm256_permute2x128_si256(
      _mm256_permute4x64_epi64((__m256i)a,216),
      _mm256_permute4x64_epi64((__m256i)b,141),48);
  }
};
template<>
struct SubsampledConcatAndCut<double,__m256d,1> {
  static __m256d  Concat( __m256d a, __m256d b, __m256d c) {
    return (__m256d) _mm256_permute2x128_si256(
      _mm256_permute4x64_epi64((__m256i)a,141),
      _mm256_permute4x64_epi64((__m256i)b,216),48);
  }
};
template<>
struct SubsampledConcatAndCut<double,__m256d,2> {
  static __m256d  Concat( __m256d a, __m256d b, __m256d c) {
    auto x = _mm256_permute2x128_si256((__m256i)a,(__m256i)c,97);
    return _mm256_blend_pd((__m256d)_mm256_permute4x64_epi64(x,180),
      (__m256d)_mm256_permute4x64_epi64((__m256i)b,225),0b0110);
  }
};
template<>
struct SubsampledConcatAndCut<double,__m256d,3> {
  static __m256d  Concat( __m256d a, __m256d b, __m256d c) {
    auto x = _mm256_permute2x128_si256((__m256i)a,(__m256i)c,97);
    return _mm256_blend_pd((__m256d)_mm256_permute4x64_epi64(x,225),
      (__m256d)_mm256_permute4x64_epi64((__m256i)b,180),0b0110);
  }
};
#elif defined USE_NEON

#endif
#endif //SUBSAMPLEDCONCATANDCUT_H
