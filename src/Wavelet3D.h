#ifndef WAVELET3D_H
#define WAVELET3D_H

// Local
#include "Wavelet.h"

// STL
#include <vector>

// Local
#include "coeffContainer.h"
#include "filters.h"
#include "vectorization.h"


/** \class Wavelet3D
 * \brief Inheritance of Wavelet class for the 1 dimensional case
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class CoeffContainerT, class WaveletSchemeT>
class Wavelet3D : public Wavelet<T,CoeffContainerT, WaveletSchemeT> {
 public:
  /// Constructor with zero initialization
  Wavelet3D() {
  
  }
  /// Constructor : Wavelet from image
  Wavelet3D(T* img, int Nc, int Nr, int Ns, bool doCycleSpinning,
      const std::string& wname, int level) : Wavelet<T,CoeffContainerT,
      WaveletSchemeT>(img, Nc, Nr, Ns, doCycleSpinning, wname, level) {

  }
  /// Default destructor
  virtual ~Wavelet3D()=default;

  /// Forward wavelet tranform
  virtual int forward() {
    return 1;
  }
  /// Backward wavelet transform: transpose of the forward transpose
  virtual int backward() {
    return 1;
  }
  /// Inverse of the wavelet tranform
  virtual int inverse() {
    return 1;
  }
};

/// Convenient type alias
template<typename T>
using PackedContainer3D =
  CoeffContainer3D<T,std::vector<T,PackAllocator<T>>>;

// Aliasing ugly types into more simple ones
template<typename T>
using Daub2_3D = Wavelet3D<T,PackedContainer3D<T>,Daub2<T>>;

/** \struct DB3DWt
 * \brief Utility struct that allow to instanciate all 3D wavelets at once
 *
 * \author Thibault Notargiacomo
 */
template<typename T>
struct DB3DWt {
 Daub2_3D<T> daub2_3D;
};



#endif //WAVELET3D_H
