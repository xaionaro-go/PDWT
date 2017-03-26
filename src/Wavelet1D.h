#ifndef WAVELET1D_H
#define WAVELET1D_H

// Local
#include "Wavelet.h"

// STL
#include <vector>

// Local
#include "coeffContainer.h"
#include "vectorization.h"

template<typename T>
using PackedContainer1D =
  CoeffContainer1D<T,std::vector<T,PackAllocator<T>>>;

/** \class Wavelet1D
 * \brief Inheritance of Wavelet class for the 1 dimensional case
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class CoeffContainerT, class WaveletSchemeT>
class Wavelet1D : Wavelet<T,CoeffContainerT, WaveletSchemeT> {
 public:
  /// Constructor with zero initialization
  Wavelet1D();
  /// Constructor : Wavelet from image
  Wavelet1D(T* img, int Nc, int Nr, int Ns, bool doCycleSpinning,
    const std::string& wname, int level);
  /// Copy constructor
  Wavelet1D(const Wavelet1D& w)=delete;
  /// Default destructor
  virtual ~Wavelet1D()=default;

  /// Forward wavelet tranform
  virtual int forward();
  /// Backward wavelet transform: transpose of the forward transpose
  virtual int backward();
  /// Inverse of the wavelet tranform
  virtual int inverse();
};

#endif //WAVELET1D_H
