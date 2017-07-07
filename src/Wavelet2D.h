#ifndef WAVELET2D_H
#define WAVELET2D_H

// Local
#include "Wavelet.h"

// STL
#include <vector>

// Local
#include "coeffContainer.h"
#include "filters.h"
#include "vectorization/vectorization.h"


/** \class Wavelet2D
 * \brief Inheritance of Wavelet class for the 1 dimensional case
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class CoeffContainerT, class WaveletSchemeT>
class Wavelet2D : public Wavelet<T,CoeffContainerT, WaveletSchemeT> {
 public:
  /// Constructor with zero initialization
  Wavelet2D()=default;

  /// Constructor : Wavelet from image
  Wavelet2D(T* img, int Nc, int Nr, int Ns, bool doCycleSpinning,
      const std::string& wname, int level) : Wavelet<T,CoeffContainerT,
      WaveletSchemeT>(img, Nc, Nr, Ns, doCycleSpinning, wname, level) {
    size_t size = Nc*Nr*Ns;
    this->m_coeff=std::make_unique<CoeffContainerT>(
      std::vector<size_t>{size}, level);
  }
  /// Default destructor
  virtual ~Wavelet2D()=default;

  /// Forward wavelet tranform
  virtual int forward() {
    /*
     * At first step, input is simply input image
     * We chose to filter the most memory friendly direction at the end in
     * order to mimic what may be done in overcomplete wavelet systems, where
     * data size grows along filtering steps
     */
/*    T* inlow = this->m_image;
    T* outlow = this->m_coeff->GetTmpBuffPtr(0,0);
    T* inhigh = this->
    for (int l=0; l<this->m_level; l++) {
      std::cout<<"Size is "<<this->m_coeff->GetScaleShape(l).at(0)<<std::endl;
      std::cout<<"Writing High subspace pointer at scale "<<l<<std::endl;
      //Y filtering
      SeparableSubsampledConvolutionEngine<T,
          typename WaveletSchemeT::f_l,
          typename WaveletSchemeT::f_h
          >::PerformSubsampledFilteringXRef(
        inlow,
        this->m_coeff->GetScaleShape(l).at(0),
        Accumulator<T,T>(outlow),
        Accumulator<T,T>(this->m_coeff->GetHighSubspacePtr(l,0)));
      //Update lowpass input and output, order is important here
      if (this->m_level==1) {
        //nothing to do
      } else if (l==this->m_level-2) {
        inlow=outlow;
        outlow=this->m_coeff->GetLowSubspacePtr(l+1);
        std::cout<<"Next scale will write to low subspace"<<std::endl;
      } else if (l==0) {
        inlow=outlow;
        outlow=this->m_coeff->GetTmpBuffPtr(l+1,0);
      } else {
        std::swap(inlow,outlow);
      }
    }*/


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
