#ifndef WAVELET2D_H
#define WAVELET2D_H

// Local
#include "Wavelet.h"

// STL
#include <vector>

// Local
#include "coeffContainer.h"
#include "filters.h"
#include "vectorization.h"


/** \class Wavelet2D
 * \brief Inheritance of Wavelet class for the 1 dimensional case
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class CoeffContainerT, class WaveletSchemeT>
class Wavelet2D : public Wavelet<T,CoeffContainerT, WaveletSchemeT> {
 public:
  /// Constructor with zero initialization
  Wavelet2D() {
  
  }
  /// Constructor : Wavelet from image
  Wavelet2D(T* img, int Nc, int Nr, int Ns, bool doCycleSpinning,
    const std::string& wname, int level) {

  }
  /// Default destructor
  virtual ~Wavelet2D()=default;

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
using PackedContainer2D =
  CoeffContainer2D<T,std::vector<T,PackAllocator<T>>>;

// Aliasing ugly types into more simple ones
template<typename T>
using Daub2_2D = Wavelet2D<T,PackedContainer2D<T>,Daub2<T>>;

/** \struct DB2DWt
 * \brief Utility struct that allow to instanciate all 2D wavelets at once
 *
 * \author Thibault Notargiacomo
 */
template<typename T>
struct DB2DWt {
 Daub2_2D<T> daub2_2D;
};

#endif //WAVELET2D_H
