#ifndef VECTORIZATION_H
#define VECTORIZATION_H

// STL
#include <vector>

// Boost
#include <boost/align/aligned_allocator.hpp>

// Local
#include "metaHelper.h"

/*
 * Include for x86 intrinsics
 * Documentation for the various intrinsics can be found on
 * https://software.intel.com/sites/landingpage/IntrinsicsGuide/
 * https://gcc.gnu.org/onlinedocs/gcc-4.8.5/gcc/ARM-NEON-Intrinsics.html
 */
#ifdef USE_SSE
  //compile using gcc -msse2
  #include "xmmintrin.h"
  #include "emmintrin.h"
#elif defined USE_AVX
  //compile using gcc -mavx2
  #include "immintrin.h"
#elif defined USE_NEON
  //compile using g++-arm-linux-gnu.x86_64 -mfpu=neon
  #include <arm_neon.h>
#endif

//Default implementation work for non-vectorized case
template<typename T, class VecT>
class VectorizedMemOp {
 public:
  static VecT load( const T* ptr ) {
    return *ptr;
  }
  static void store( float* ptr, VecT value) {
    *ptr = value;
  }
};

#ifdef USE_SSE
template<>
class VectorizedMemOp<float,__m128> {
 public:
  static __m128 load( const float* ptr ) {
    return _mm_load_ps( ptr );
  }
  static void store( float* ptr, __m128 value) {
    _mm_store_ps( ptr, value );
  }
};
#elif defined USE_AVX
template<>
class VectorizedMemOp<float,__m256> {
 public:
  static __m256 load( const float* ptr ) {
    return _mm256_load_ps( ptr );
  }
  static void store( float* ptr, __m256 value) {
    _mm256_store_ps( ptr, value );
  }
};
#elif defined USE_NEON
template<>
class VectorizedMemOp<float,float32x4_t> {
 public:
  static float32x4_t load( const float* ptr ) {
    return vld1q_f32( ptr );
  }
  static void store( float* ptr, float32x4_t value) {
    vst1q_f32( ptr, value );
  }
};
#endif

///Perform left and right shift
template<typename T, class VecT, unsigned char SHIFT>
class VectorizedShift {
 public:
  //Defaulted implementation for scalar type: no shift
  constexpr static VecT LeftShift( VecT input ) {
    return SHIFT != 0 ? 0 : input;
  }
  //Defaulted implementation for scalar type: no shift
  constexpr static VecT RightShift( VecT input ) {
    return SHIFT != 0 ? 0 : input;
  }
};

/*
 * Concatenate two vectors that are given as input
 * then right shift the results of RIGHT_SHIFT elements
 */
template<typename T, class VecT, int RIGHT_SHIFT>
class VectorizedConcatAndCut {
 public:
  //Defaulted implementation is to rely on individual shifting
  //and add
  static VecT Concat( VecT left, VecT right ) {
    //Perform shift on both operand
    //The left one should be right shifted and right one should be left shifted
    VecT r = VectorizedShift<T,VecT,RightShiftBytes>::RightShift(left);
    VecT l = VectorizedShift<T,VecT,LeftShiftBytes>::LeftShift(right);
    //return sum of the two
    return r + l;
}
 protected:
  constexpr static int VecSize = sizeof(VecT)/sizeof(T);
  constexpr static int RightShiftBytes = RIGHT_SHIFT*sizeof(T);
  constexpr static int LeftShiftBytes = (VecSize-RIGHT_SHIFT)*sizeof(T);
};

#ifdef USE_SSE
template<unsigned char SHIFT>
class VectorizedShift<float,__m128,SHIFT> {
 public:
  constexpr static __m128 LeftShift( __m128 input ) {
    return (__m128)_mm_slli_si128( (__m128i)input, SHIFT );
  }
  constexpr static __m128 RightShift( __m128 input ) {
    return (__m128)_mm_srli_si128( (__m128i)input, SHIFT );
  }
};
#elif defined USE_AVX
template<int Val, class enable=void>
struct AVX256ConcatandCut {
  static __m256 Concat(__m256 left, __m256 right) {
    assert(("Vectorized Shift AVX256 cannot account for shift > 256 bits",
          false));
    return left;
  }
};

template<int Val>
struct AVX256ConcatandCut<Val, typename ctrange<0, 1, Val>::enabled> {
  static __m256 Concat(__m256 left, __m256 right) {
    return left;
  }
};
template<int Val>
struct AVX256ConcatandCut<Val, typename ctrange<1, 4, Val>::enabled> {
  static __m256 Concat(__m256 left, __m256 right) {
    right=(__m256)_mm256_permute2x128_si256((__m256i)left,(__m256i)right,33);
    return (__m256)_mm256_alignr_epi8((__m256i)right,(__m256i)left,
        Val*sizeof(int));
  }
};
template<int Val>
struct AVX256ConcatandCut<Val, typename ctrange<4, 5, Val>::enabled> {
  static __m256 Concat(__m256 left, __m256 right) {
    return (__m256)_mm256_permute2x128_si256((__m256i)left,(__m256i)right,33);
  }
};
template<int Val>
struct AVX256ConcatandCut<Val, typename ctrange<5, 8, Val>::enabled> {
  static __m256 Concat(__m256 left, __m256 right) {
    left=(__m256)_mm256_permute2x128_si256((__m256i)left,(__m256i)right,33);
    return (__m256)_mm256_alignr_epi8((__m256i)right,(__m256i)left,
        (Val-4)*sizeof(int));
  }
};

template<int Val>
struct AVX256ConcatandCut<Val, typename ctrange<8, 9, Val>::enabled> {
  __m256 static Concat(__m256 left, __m256 right) {
    return right;
  }
};

template<int RIGHT_SHIFT>
class VectorizedConcatAndCut<float,__m256,RIGHT_SHIFT> {
 public:
  //Optimized specific intrinsic for concat / shift / cut in AVX
  static __m256 Concat( __m256 left, __m256 right ) {
    return AVX256ConcatandCut<RIGHT_SHIFT>::Concat(left,right);
  }
};
#elif defined USE_NEON
template<int RIGHT_SHIFT>
class VectorizedConcatAndCut<float,float32x4_t,RIGHT_SHIFT> {
 public:
  //Optimized specific intrinsic for concat / shift / cut in AVX
  static float32x4_t Concat( float32x4_t left, float32x4_t right ) {
    return vextq_f32( left, right, RIGHT_SHIFT) ;
  }
};
#endif

//Specialize Packed types when they exist
//Default packed type is... not packed
template<typename T> struct PackedType { using type=T; };
#ifdef USE_SSE
  template<> struct PackedType<float> { using type = __m128; };
#elif defined USE_AVX
  template<> struct PackedType<float> { using type = __m256; };
#elif defined USE_NEON
  template<> struct PackedType<float> { using type = float32x4_t; };
#endif

template<typename T> using PackType = typename PackedType<T>::type;

template<typename T>
using PackAllocator =
  boost::alignment::aligned_allocator<T,sizeof(PackType<T>)>;

#endif /* VECTORIZATION_H */
