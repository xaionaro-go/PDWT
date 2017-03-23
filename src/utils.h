#ifndef UTILS_H
#define UTILS_H

// STL
#include <type_traits>


/** \struct w_info
 * \brief Description of the workload
 *
 * \author Pierre Paleo
 */
struct w_info {
  /// Number of dimensions. For now only 2D and (batched) 1D are supported
  int ndims;
  /// Number of rows of the image (for 1D : Nr = 1)
  int Nr;
  /// Number of rows of the image (for 1D : Nr = 1)
  int Nc;
  /// Number of decomposition levels
  int nlevels;
  /// Do Stationary (Undecimated) Wavelet Transform
  int do_swt;              
  /// "Filter" length
  int hlen;
};

template<typename T, typename U>
constexpr auto w_iDivUp(T a, U b) {
  static_assert(std::is_integral<T>::value, "Integer required.");
  static_assert(std::is_integral<U>::value, "Integer required.");
  return (a % b != 0) ? (a / b + 1) : (a / b);
}

template<typename T>
constexpr auto w_ipow2(T a) {
  static_assert(std::is_integral<T>::value, "Integer required.");
  return 1 << a;
}

template<typename T>
constexpr auto w_ilog2(T i) {
  static_assert(std::is_integral<T>::value, "Integer required.");
  int l = 0;
  while (i >>= 1) {
    ++l;
  }
  return l;
}


/// When the size is odd, allocate one extra element before subsampling
template<typename T>
void w_div2(T* N) {
  static_assert(std::is_integral<T>::value, "Integer required.");
  if ((*N) & 1)
    *N = ((*N)+1)/2;
  else
    *N = (*N)/2;
}

#endif
