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
  static int PerformSubsampledFilteringC(const T* in,
      T* out, w_info info) {
    
    int Nc=w_info.Nc; 
    int Nc_is_odd = (Nc & 1);
    int NxOut = (Nc + Nc_is_odd)/2;

    // Loop over output image
    for (int ox=0; ox<NxOut; ox++) {
      T acc=0;
      // Loop over filter size, with periodic boundary conditions
      for (int fx=-FilterT::TapSizeLeft; fx<FilterT::TapSizeRight; fx++) {
        int ix = ox*2 + fx;
        // if N is odd, image is virtually extended
        if (ix < 0) ix += (Nc + Nc_is_odd);
        // no "else if", since idx_x can be > N-1  after being incremented
        if (ix > Nc-1) {
          // if N is odd, repeat the right-most element
          if ((ix == Nc) && (Nc_is_odd))
            ix--;
          // if N is odd, image is virtually extended
          else
            ix -= (Nc + Nc_is_odd);
        }
        acc += in[ix] * FilterT::Buf[fx];
      }
      out[ox] = acc;
    }
  }
};

#endif //SEPARABLE_H
