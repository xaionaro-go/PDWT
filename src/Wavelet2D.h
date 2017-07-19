#ifndef WAVELET2D_H
#define WAVELET2D_H

// Local
#include "Wavelet.h"

// STL
#include <vector>

// Local
#include "coeffContainer.h"
#include "filters.h"
#include "separable.h"
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
    /**
     * At first step, input is simply input image
     * We chose to filter the most memory friendly direction at the end in
     * order to mimic what may be done in overcomplete wavelet systems, where
     * data size grows along filtering steps
     */
    T* inlow = this->m_image;
    T* outlow = this->m_coeff->GetHalfTmpBuffPtr(0);
    T* outhigh = this->m_coeff->GetHalfTmpBuffPtr(1);
    for (int l=0; l<this->m_level; l++) {
      //Y filtering
      SeparableSubsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::f_l,
          typename WaveletSchemeT::f_h
          >::PerformSubsampledFilteringYRef(
        this->m_coeff->GetScaleShape(l).at(0),
        this->m_coeff->GetScaleShape(l).at(1),
        Accumulator<T,T,T,int>(inlow, outlow,
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l+1).at(0)),
        Accumulator<T,T,T,int>(inlow, outhigh,
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l+1).at(0)));

      T* inlowY = outlow;
      T* inhighY = outhigh;

      if (this->m_level==1) {
        outlow=this->m_coeff->GetLowSubspacePtr(l);
      } else {
        outlow=this->m_coeff->GetOutLowTmpBuffPtr();
      }

      //Now perform X filtering on lowpass Y
      SeparableSubsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::f_l,
          typename WaveletSchemeT::f_h
          >::PerformSubsampledFilteringXRef(
        this->m_coeff->GetScaleShape(l).at(0),
        this->m_coeff->GetScaleShape(l+1).at(1),
        Accumulator<T,T,T,int>(inlowY, outlow),
        Accumulator<T,T,T,int>(inlowY,
          this->m_coeff->GetHighSubspacePtr(l,1)));
      //Now perform X filtering on highpass Y
      SeparableSubsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::f_l,
          typename WaveletSchemeT::f_h
          >::PerformSubsampledFilteringXRef(
        this->m_coeff->GetScaleShape(l).at(0),
        this->m_coeff->GetScaleShape(l+1).at(1),
        Accumulator<T,T,T,int>(inhighY,
          this->m_coeff->GetHighSubspacePtr(l,2)),
        Accumulator<T,T,T,int>(inhighY,
          this->m_coeff->GetHighSubspacePtr(l,3)));

      //Update lowpass input and output, order is important here
      inlow=outlow;
    }
    return 1;
  }
  /// Backward wavelet transform: transpose of the forward transform
  virtual int backward() {
    for (int l=this->m_level; l>0; l--) {
      // Invert X lowpass/highpass filtering for lowpass Y
/*      SeparableUpsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::i_l,
          typename WaveletSchemeT::i_h
        >::PerformUpsampledFilteringXRef(
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l-1).at(0),
          this->m_coeff->GetScaleShape(l).at(1),
          this->m_coeff->GetScaleShape(l-1).at(1),
          SubsampledAccumulator<T,T,int,int,
              typename WaveletSchemeT::i_l,
              typename WaveletSchemeT::i_h>(
            this->m_coeff->GetHalfTmpBuffPtr(0)),
          this->m_coeff->GetLowSubspacePtr(l-1),
          this->m_coeff->GetHighSubspacePtr(l-1,1));

      // Invert X lowpass/highpass filtering for highpass Y
      SeparableUpsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::i_l,
          typename WaveletSchemeT::i_h
        >::PerformUpsampledFilteringXRef(
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l-1).at(0),
          this->m_coeff->GetScaleShape(l).at(1),
          this->m_coeff->GetScaleShape(l-1).at(1),
          SubsampledAccumulator<T,T,int,int,
              typename WaveletSchemeT::i_l,
              typename WaveletSchemeT::i_h>(
          this->m_coeff->GetHalfTmpBuffPtr(1)),
          this->m_coeff->GetHighSubspacePtr(l-1,2),
          this->m_coeff->GetHighSubspacePtr(l-1,3));

      //Update output buffer destination
      T* outlow;
      if (l<=1) {
        outlow=this->m_image;
      } else {
        outlow=this->m_coeff->GetLowSubspacePtr(l);
      }

      // Invert Y lowpass/highpass filtering
      SeparableUpsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::i_l,
          typename WaveletSchemeT::i_h
        >::PerformUpsampledFilteringYRef(
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l-1).at(0),
          this->m_coeff->GetScaleShape(l).at(1),
          this->m_coeff->GetScaleShape(l-1).at(1),
          SubsampledAccumulator<T,T,int,int,
              typename WaveletSchemeT::i_l,
              typename WaveletSchemeT::i_h>(outlow),
            this->m_coeff->GetHalfTmpBuffPtr(0),
            this->m_coeff->GetHalfTmpBuffPtr(1));*/
    }
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
template<typename T>
using Daub3_2D = Wavelet2D<T,PackedContainer2D<T>,Daub3<T>>;
template<typename T>
using Daub4_2D = Wavelet2D<T,PackedContainer2D<T>,Daub4<T>>;
template<typename T>
using Daub5_2D = Wavelet2D<T,PackedContainer2D<T>,Daub5<T>>;
template<typename T>
using Anto97_BiOrth_2D = Wavelet2D<T,PackedContainer2D<T>,Anto97_BiOrth<T>>;
template<typename T>
using QSHIFT6_Orth_2D = Wavelet2D<T,PackedContainer2D<T>,QSHIFT6_Orth<T>>;
template<typename T>
using REVERSE_QSHIFT6_Orth_2D = 
  Wavelet2D<T,PackedContainer2D<T>,REVERSE_QSHIFT6_Orth<T>>;


// Aliasing ugly types into more simple ones
//template<typename T>
//using dtwAnto97QSHIFT6_2D = 
//  DTWavelet2D<T,PackedDTContainer2D<T>,dtwAnto97QSHIFT6<T>>;

/** \struct DB2DWt
 * \brief Utility struct that allow to instanciate all 2D wavelets at once
 *
 * \author Thibault Notargiacomo
 */
template<typename T>
struct DB2DWt {
 Daub2_2D<T> daub2_2D;
 Daub3_2D<T> daub3_2D;
 Daub4_2D<T> daub4_2D;
 Daub5_2D<T> daub5_2D;
 Anto97_BiOrth_2D<T> anto97_BiOrth_2D;
 QSHIFT6_Orth_2D<T> QShift6_Orth_2D;
 REVERSE_QSHIFT6_Orth_2D<T> Reverse_Qshift6_Orth_2D;
// dtwAnto97QSHIFT6_2D<T> dtwAnto97QShift6_2D; 
};

#endif //WAVELET2D_H
