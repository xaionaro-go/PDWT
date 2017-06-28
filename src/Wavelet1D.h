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
#include "vectorization/vectorization.h"


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
    size_t size = Nc*Nr*Ns;
    this->m_coeff=std::make_unique<CoeffContainerT>(
      std::vector<size_t>{size}, level);
  }
  /// Default destructor
  virtual ~Wavelet1D()=default;

  /// Forward wavelet tranform
  virtual int forward() {
    // At first step, input is simply input image
    T* inlow = this->m_image;
    // Lowpass output is either in the wv tree if we have already finished,
    // or it is stored to a temporary buffer for further decomposition
    T* outlow = (this->m_level<2) ? this->m_coeff->GetLowSubspacePtr(1) :
      this->m_coeff->GetTmpBuffPtr().at(0)->data();
    for (int l=0; l<this->m_level; l++) {
      std::cout<<"Size is "<<this->m_coeff->GetScaleShape(l).at(0)<<std::endl;
      std::cout<<"Writing High subspace pointer at scale "<<l<<std::endl;
      //#pragma omp parallel for
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
        outlow=this->m_coeff->GetTmpBuffPtr().at(1)->data();
      } else {
        std::swap(inlow,outlow);
      }
    }
    return 1;
  }
  /// Backward wavelet transform: transpose of the forward transform
  virtual int backward() {
    // At first step, input lowpass image is part of the wavelet tree
    T* inlow = this->m_coeff->GetLowSubspacePtr(this->m_level-1);
    // output is either in the image if we have already finished,
    // or if number of level is odd. Otherwise, it is stored to a
    //  temporary buffer for further reconstruction
    T* outlow = ((this->m_level%2==1) ? this->m_image :
      this->m_coeff->GetTmpBuffPtr().at(0)->data());
    for (int l=this->m_level; l>0; l--) {
      std::cout<<"Inverse, taking as input 2 subspaces of size "<<
        this->m_coeff->GetScaleShape(l).at(0)<<" In order to form a new"<<
        " image of size "<<this->m_coeff->GetScaleShape(l-1).at(0)<<std::endl;
      //std::cout<<"BufSize is "<<this->m_coeff->GetTmpBuffPtr().at(0)->size()
      //  <<std::endl;
      std::cout<<"Current scale is: at scale "<<l<<std::endl;
      //#pragma omp parallel for
      SeparableUpsampledConvolutionEngine<T,
          typename WaveletSchemeT::i_l,
          typename WaveletSchemeT::i_h
        >::PerformUpsampledFilteringXRef(
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l-1).at(0),
          outlow,
          inlow,
          this->m_coeff->GetHighSubspacePtr(l-1,0));

      //Update lowpass input and output
      if (l==this->m_level && l>1) {
        inlow= ((this->m_level%2==0) ? this->m_image :
          this->m_coeff->GetTmpBuffPtr().at(0)->data());
       }
       std::swap(inlow,outlow);
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
using PackedContainer1D =
  CoeffContainer1D<T,std::vector<T>>;

// Aliasing ugly types into more simple ones
template<typename T>
using Daub2_1D = Wavelet1D<T,PackedContainer1D<T>,Daub2<T>>;
template<typename T>
using Daub3_1D = Wavelet1D<T,PackedContainer1D<T>,Daub3<T>>;
template<typename T>
using Daub4_1D = Wavelet1D<T,PackedContainer1D<T>,Daub4<T>>;
template<typename T>
using Daub5_1D = Wavelet1D<T,PackedContainer1D<T>,Daub5<T>>;
template<typename T>
using Anto97_BiOrth_1D = Wavelet1D<T,PackedContainer1D<T>,Anto97_BiOrth<T>>;
template<typename T>
using QSHIFT6_Orth_1D = Wavelet1D<T,PackedContainer1D<T>,QSHIFT6_Orth<T>>;
template<typename T>
using REVERSE_QSHIFT6_Orth_1D = 
  Wavelet1D<T,PackedContainer1D<T>,REVERSE_QSHIFT6_Orth<T>>;


/** \class DTWavelet1D
 * \brief Inheritance of Wavelet class for the Dual Tree 1 dimensional case
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class CoeffContainerT, class DTWaveletSchemeT>
class DTWavelet1D : public Wavelet1D<T,CoeffContainerT, DTWaveletSchemeT> {
 public:
  /// Constructor with zero initialization
  DTWavelet1D()=default;
  
  /// Constructor : Wavelet from image
  DTWavelet1D(T* img, int Nc, int Nr, int Ns, bool doCycleSpinning,
      const std::string& wname, int level) : Wavelet1D<T,CoeffContainerT,
      DTWaveletSchemeT>(img, Nc, Nr, Ns, doCycleSpinning, wname, level) {}

  /// Default destructor
  virtual ~DTWavelet1D()=default;

  /// Forward wavelet tranform
  virtual int forward() {
    // At first step, input is simply input image
    T* inlow = this->m_image;
    // Lowpass output is either in the wv tree if we have already finished,
    // or it is stored to a temporary buffer for further decomposition
    T* outlowReal = (this->m_level<2) ? this->m_coeff->GetLowSubspacePtr(1,0) :
      this->m_coeff->GetTmpBuffPtr().at(0)->data();
    T* outlowImag = (this->m_level<2) ? this->m_coeff->GetLowSubspacePtr(1,1) :
      this->m_coeff->GetTmpBuffPtr().at(0)->data();


    for (int l=0; l<this->m_level; l++) {
      std::cout<<"Size is "<<this->m_coeff->GetScaleShape(l).at(0)<<std::endl;
      std::cout<<"Writing High subspace pointer at scale "<<l<<std::endl;
      //#pragma omp parallel for
      SeparableSubsampledConvolutionEngine<T,
          typename DTWaveletSchemeT::f_l0r,
          typename DTWaveletSchemeT::f_h0r
          >::PerformSubsampledFilteringXRef(
        inlow,
        this->m_coeff->GetScaleShape(l).at(0),
        outlowReal,
        this->m_coeff->GetHighSubspacePtr(l,0));
      //Update lowpass input and output, order is important here
      if (this->m_level==1) {
        //nothing to do
      } else if (l==this->m_level-2) {
        inlow=outlowReal;
        outlowReal=this->m_coeff->GetLowSubspacePtr(l+1);
        std::cout<<"Next scale will write to low subspace"<<std::endl;
      } else if (l==0) {
        inlow=outlowReal;
        outlowReal=this->m_coeff->GetTmpBuffPtr().at(1)->data();
      } else {
        std::swap(inlow,outlowReal);
      }
    }
    return 1;
  }
  /// Backward wavelet transform: transpose of the forward transform
  virtual int backward() {
    // At first step, input lowpass image is part of the wavelet tree
    T* inlow = this->m_coeff->GetLowSubspacePtr(this->m_level-1);
    // output is either in the image if we have already finished,
    // or if number of level is odd. Otherwise, it is stored to a
    //  temporary buffer for further reconstruction
    T* outlow = ((this->m_level%2==1) ? this->m_image :
      this->m_coeff->GetTmpBuffPtr().at(0)->data());
    for (int l=this->m_level; l>0; l--) {
      std::cout<<"Inverse, taking as input 2 subspaces of size "<<
        this->m_coeff->GetScaleShape(l).at(0)<<" In order to form a new"<<
        " image of size "<<this->m_coeff->GetScaleShape(l-1).at(0)<<std::endl;
      //std::cout<<"BufSize is "<<this->m_coeff->GetTmpBuffPtr().at(0)->size()
      //  <<std::endl;
      std::cout<<"Current scale is: at scale "<<l<<std::endl;
      //#pragma omp parallel for
      SeparableUpsampledConvolutionEngine<T,
          typename DTWaveletSchemeT::i_l0r,
          typename DTWaveletSchemeT::i_h0r
        >::PerformUpsampledFilteringXRef(
          inlow,
          this->m_coeff->GetHighSubspacePtr(l-1,0),
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l-1).at(0),
          outlow);

      //Update lowpass input and output
      if (l==this->m_level && l>1) {
        inlow= ((this->m_level%2==0) ? this->m_image :
          this->m_coeff->GetTmpBuffPtr().at(0)->data());
       }
       std::swap(inlow,outlow);
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
using PackedDTContainer1D =
  DTCoeffContainer1D<T,std::vector<T>>;

// Aliasing ugly types into more simple ones
template<typename T>
using dtwAnto97QSHIFT6_1D = 
  DTWavelet1D<T,PackedDTContainer1D<T>,dtwAnto97QSHIFT6<T>>;

/** \struct All1DWavelet
 * \brief Utility struct that allow to instanciate all 1D wavelets at once
 *
 * \author Thibault Notargiacomo
 */
template<typename T>
struct DB1DWt {
 Daub2_1D<T> daub2_1D;
 Daub3_1D<T> daub3_1D;
 Daub4_1D<T> daub4_1D;
 Daub5_1D<T> daub5_1D;
 Anto97_BiOrth_1D<T> anto97_BiOrth_1D;
 QSHIFT6_Orth_1D<T> QShift6_Orth_1D;
 REVERSE_QSHIFT6_Orth_1D<T> Reverse_Qshift6_Orth_1D;
 dtwAnto97QSHIFT6_1D<T> dtwAnto97QShift6_1D; 
};

/*
-Cas 1 seul niveau :
  LP et HP on ete ecrits direct dans les coeff
  la reconstruction peut se faire direct dans image

-Cas 2 niveaux :
  On a alloue un buffer temp0 de taille lev1 pour mettre le LP1
  coeff contient HP1 HP2 LP2
  Lors de la reco on a besoin de reconstruire LP1 dans temp0

-Cas 3 niveaux :
  On alloue un buffer temp0 de taille lev1 pour mettre le LP1
  coeff contient HP1 HP2 HP3 LP3
  Lors de la reco on a besoin de reconstruire LP2 dans image
  Puis on reco LP1 dans temp0

-Cas 4 niveaux:
  On alloue un buffer temp0 de taille lev1 pour mettre le LP1
  coeff contient HP1 HP2 HP3 HP4 LP4
  Lors de la reco on a besoin de reconstruire LP3 dans temp0
  Puis on reco LP2 dans image
  Puis on reco LP1 dans temp0
*/
#endif //WAVELET1D_H
