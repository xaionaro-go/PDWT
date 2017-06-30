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
      this->m_coeff->GetTmpBuffPtr(0,0);
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
        outlow=this->m_coeff->GetTmpBuffPtr(l+1,0);
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
    // In 1D, one has to keep in mind, that level 1 (last recontruction step)
    // the input low image has to be of size image/2.
    // this means that, at start, odd level must use the biggest tmp buff
    T* outlow =  this->m_level > 1 ? 
      this->m_coeff->GetTmpBuffPtr(this->m_level,0) : this->m_image;

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
      if (l<=2) {
        inlow = outlow;
        outlow = this->m_image;
       } else if (l==this->m_level) {
         inlow = outlow;
         outlow = this->m_coeff->GetTmpBuffPtr(this->m_level-1,0);
       } else {
         std::swap(inlow,outlow);
       }
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
class DTWavelet1D : public Wavelet<T,CoeffContainerT, DTWaveletSchemeT> {
 public:
  /// Constructor with zero initialization
  DTWavelet1D()=default;
  
  /// Constructor : Wavelet from image
  DTWavelet1D(T* img, int Nc, int Nr, int Ns, bool doCycleSpinning,
      const std::string& wname, int level) : Wavelet<T,CoeffContainerT,
      DTWaveletSchemeT>(img, Nc, Nr, Ns, doCycleSpinning, wname, level) {
    size_t size = Nc*Nr*Ns;
    this->m_coeff=std::make_unique<CoeffContainerT>(
      std::vector<size_t>{size}, level);
  }

  /// Default destructor
  virtual ~DTWavelet1D()=default;

  /// Forward wavelet tranform
  virtual int forward() {
    // At first step, input is simply input image
    T* inlow = this->m_image;
    T* inlowReal = nullptr;
    T* inlowImag = nullptr;

    // Lowpass output is either in the wv tree if we have already finished,
    // or it is stored to a temporary buffer for further decomposition
    T* outlowReal = (this->m_level<2) ? this->m_coeff->GetLowSubspacePtr(1,0) :
      this->m_coeff->GetTmpBuffPtr(0,0);
    //output size is size/2, so one has to get proper offset
    T* outlowImag = (this->m_level<2) ? this->m_coeff->GetLowSubspacePtr(1,1) :
      this->m_coeff->GetTmpBuffPtr(0,1);

    // Divide the tree in two branch along the X direction
    int l=0;
    if (this->m_level>=1) {
      //STRONG TODO: use is_equal metaprogramming technique to check if filters
      // are the same, and perform much less work !!!!
      SeparableSubsampledConvolutionEngine<T,
          typename DTWaveletSchemeT::f_l0r,
          typename DTWaveletSchemeT::f_l0i,
          typename DTWaveletSchemeT::f_h0r,
          typename DTWaveletSchemeT::f_h0i
          >::PerformSubsampledFilteringXRef(
        inlow,
        this->m_coeff->GetScaleShape(l).at(0),
        Accumulator<T,T>(outlowReal),
        Accumulator<T,T>(outlowImag),
        Accumulator<T,T>(this->m_coeff->GetHighSubspacePtr(l,0,0)),
        Accumulator<T,T>(this->m_coeff->GetHighSubspacePtr(l,0,1)));

      //Update lowpass input and output, order is important here
      if (this->m_level==2) {
        inlowReal=outlowReal;
        inlowImag=outlowImag;
        outlowReal=this->m_coeff->GetLowSubspacePtr(this->m_level,0);
        outlowImag=this->m_coeff->GetLowSubspacePtr(this->m_level,1);
      } else if (this->m_level>2) {
        inlowReal=outlowReal;
        inlowImag=outlowImag;
        outlowReal=this->m_coeff->GetTmpBuffPtr(1,0);
        outlowImag=this->m_coeff->GetTmpBuffPtr(1,1);
      }
    }

    //Band 0: real X
    for (int l=1; l<this->m_level; l++) {
      SeparableSubsampledConvolutionEngine<T,
          typename DTWaveletSchemeT::f_lnr,
          typename DTWaveletSchemeT::f_hnr
          >::PerformSubsampledFilteringXRef(
        inlowReal,
        this->m_coeff->GetScaleShape(l).at(0),
        Accumulator<T,T>(outlowReal),
        Accumulator<T,T>(this->m_coeff->GetHighSubspacePtr(l,0,0)));
      //Update lowpass input and output, order is important here
      if (l>=this->m_level-2) {
        inlowReal=outlowReal;
        outlowReal=this->m_coeff->GetLowSubspacePtr(this->m_level,0);
      } else {
        std::swap(inlowReal,outlowReal);
      }
    }
    //Band 1: imag X
    for (int l=1; l<this->m_level; l++) {
      SeparableSubsampledConvolutionEngine<T,
          typename DTWaveletSchemeT::f_lni,
          typename DTWaveletSchemeT::f_hni
          >::PerformSubsampledFilteringXRef(
        inlowImag,
        this->m_coeff->GetScaleShape(l).at(0),
        Accumulator<T,T>(outlowImag),
        Accumulator<T,T>(this->m_coeff->GetHighSubspacePtr(l,0,1)));
      //Update lowpass input and output, order is important here
      if (l>=this->m_level-2) {
        inlowImag=outlowImag;
        outlowImag=this->m_coeff->GetLowSubspacePtr(this->m_level,1);
      } else {
        std::swap(inlowImag,outlowImag);
      }
    }
    // map the set of filtered signals to the real DTCWT mixture
    this->m_coeff->WaveletToCpx();
    return 1;
  }

  /// Backward wavelet transform: transpose of the forward transform
  virtual int backward() {

    // map the set of DTCWT coefficients to a simple set of filtered signals
    this->m_coeff->CpxToWavelet();
 
    // Band 0: real X
    T* inlowReal = this->m_coeff->GetLowSubspacePtr(this->m_level-1,0);
    T* outlowReal= this->m_coeff->GetTmpBuffPtr(this->m_level,0);
    for (int l=this->m_level; l>1; l--) {
      //#pragma omp parallel for
      SeparableUpsampledConvolutionEngine<T,
          typename DTWaveletSchemeT::i_lnr,
          typename DTWaveletSchemeT::i_hnr
        >::PerformUpsampledFilteringXRef(
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l-1).at(0),
          outlowReal,
          inlowReal,
          this->m_coeff->GetHighSubspacePtr(l-1,0,0));

      //Update lowpass input and output
      if (l==this->m_level) {
        inlowReal=this->m_coeff->GetTmpBuffPtr(this->m_level-1,0);
      }
      std::swap(inlowReal,outlowReal);
    }
    auto print = [](auto in){std::cout<<in<<" , ";};
    std::cout<<"Print low real: ";
    std::for_each(inlowReal,inlowReal+5,print);
    std::cout<<std::endl;

   // Band 1: imag X
    T* inlowImag = this->m_coeff->GetLowSubspacePtr(this->m_level-1,1);
    T* outlowImag= this->m_coeff->GetTmpBuffPtr(this->m_level,1);
    for (int l=this->m_level; l>1; l--) {
	  //#pragma omp parallel for
	  SeparableUpsampledConvolutionEngine<T,
		  typename DTWaveletSchemeT::i_lni,
		  typename DTWaveletSchemeT::i_hni
		>::PerformUpsampledFilteringXRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l-1).at(0),
		  outlowImag,
		  inlowImag,
		  this->m_coeff->GetHighSubspacePtr(l-1,0,1));
	  //Update lowpass input and output
      if (l==this->m_level) {
        inlowImag=this->m_coeff->GetTmpBuffPtr(this->m_level-1,1);
      }
	  std::swap(inlowImag,outlowImag);
	}
    std::cout<<"Print low imag: ";
    std::for_each(inlowImag,inlowImag+5,print);
    std::cout<<std::endl;

    // Merge both trees in the X direction
    T* outlow = this->m_image;
    for (int l=this->m_level; l>0; l--) {
      //#pragma omp parallel for
      SeparableUpsampledConvolutionEngine<T,
          typename DTWaveletSchemeT::i_l0r,
          typename DTWaveletSchemeT::i_l0i,
          typename DTWaveletSchemeT::i_h0r,
          typename DTWaveletSchemeT::i_h0i
        >::PerformUpsampledFilteringXRef(
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l-1).at(0),
          outlow,
          inlowReal,
          inlowImag,
          this->m_coeff->GetHighSubspacePtr(l-1,0,0),
          this->m_coeff->GetHighSubspacePtr(l-1,0,1));
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
