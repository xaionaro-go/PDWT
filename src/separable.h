#ifndef SEPARABLE_H
#define SEPARABLE_H

// STL

// Local
#include "filters.h"

/** \class SeparableSubsampledConvolutionEngine
 * \brief Code for the separable subsample convolution
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class FilterT>
class SeparableSubsampledConvolutionEngine {
 public:
  /// Defaulted constructor
  SeparableSubsampledConvolutionEngine()=default;
  /// Default destructor
  virtual ~SeparableSubsampledConvolutionEngine()=default;

  /// The main method : perform Subsampled convolution on one row
  static int PerformSubsampledFilteringR(const T* image,
    T* out, w_info info) {
 /* int gidx = threadIdx.x + blockIdx.x*blockDim.x;
  int gidy = threadIdx.y + blockIdx.y*blockDim.y;
  int Nc_is_odd = (Nc & 1);
  int Nc2 = (Nc + Nc_is_odd)/2;

  // horiz subsampling : Input (Nr, Nc) => Output (Nr, Nc/2)
  if (gidy < Nr && gidx < Nc2) {
    int c, hL, hR;
    // odd kernel size
    if (hlen & 1) {
      c = hlen/2;
      hL = c;
      hR = c;
    }
    else { // even kernel size : center is shifted to the left
      c = hlen/2 - 1;
      hL = c;
      hR = c+1;
    }
    DTYPE res_tmp_a1 = 0, res_tmp_a2 = 0;
    DTYPE img_val;

    // Convolution with periodic boundaries extension.
    for (int jx = 0; jx <= hR+hL; jx++) {
      int idx_x = gidx*2 - c + jx;
      // if N is odd, image is virtually extended
      if (idx_x < 0) idx_x += (Nc + Nc_is_odd);
      // no "else if", since idx_x can be > N-1  after being incremented
      if (idx_x > Nc-1) {
        // if N is odd, repeat the right-most element
        if ((idx_x == Nc) && (Nc_is_odd))
          idx_x--;
        // if N is odd, image is virtually extended
        else
          idx_x -= (Nc + Nc_is_odd);
      }
      img_val = img[gidy*Nc + idx_x];
      res_tmp_a1 += img_val * c_kern_L[hlen-1 - jx];
      res_tmp_a2 += img_val * c_kern_H[hlen-1 - jx];
    }
    tmp_a1[gidy* Nc2 + gidx] = res_tmp_a1;
    tmp_a2[gidy* Nc2 + gidx] = res_tmp_a2;
 */
    return 0;
  }
};

#endif //SEPARABLE_H
