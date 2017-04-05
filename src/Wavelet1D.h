#ifndef WAVELET1D_H
#define WAVELET1D_H

// Local
#include "Wavelet.h"

// STL
#include <vector>

// Local
#include "coeffContainer.h"
#include "filters.h"
#include "separable.h"
#include "vectorization.h"


/** \class Wavelet1D
 * \brief Inheritance of Wavelet class for the 1 dimensional case
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class CoeffContainerT, class WaveletSchemeT>
class Wavelet1D : public Wavelet<T,CoeffContainerT, WaveletSchemeT> {
 public:
  /// Constructor with zero initialization
  Wavelet1D()=default;
  
  /// Constructor : Wavelet from image
  Wavelet1D(T* img, int Nc, int Nr, int Ns, bool doCycleSpinning,
      const std::string& wname, int level) : Wavelet<T,CoeffContainerT,
      WaveletSchemeT>(img, Nc, Nr, Ns, doCycleSpinning, wname, level) {
    size_t size = Nc*Nc*Nr;
    this->m_coeff=std::make_unique<CoeffContainerT>(
      std::vector<size_t>{size}, level);
  }
  /// Default destructor
  virtual ~Wavelet1D()=default;

  /// Forward wavelet tranform
  virtual int forward() {
    int Nx=this->m_info.Nx;
    T* in = this->m_image;
    for (int l=0; l<this->m_level; l++) {
      //#pragma omp parallel for
      SeparableSubsampledConvolutionEngine<T,
          typename WaveletSchemeT::f_l,
          typename WaveletSchemeT::f_h
        >::PerformSubsampledFilteringXRef(
          in, Nx,
          this->m_coeff->GetLowSubspacePtr(l),
          this->m_coeff->GetHighSubspacePtr(l,0));
       //Update lowpass input
       in=this->m_coeff->GetLowSubspacePtr(l);
    }
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
using PackedContainer1D =
  CoeffContainer1D<T,std::vector<T,PackAllocator<T>>>;

// Aliasing ugly types into more simple ones
template<typename T>
using Daub2_1D = Wavelet1D<T,PackedContainer1D<T>,Daub2<T>>;

/** \struct All1DWavelet
 * \brief Utility struct that allow to instanciate all 1D wavelets at once
 *
 * \author Thibault Notargiacomo
 */
template<typename T>
struct DB1DWt {
 Daub2_1D<T> daub2_1D;
};

#endif //WAVELET1D_H
