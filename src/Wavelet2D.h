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
    this->m_coeff=std::make_unique<CoeffContainerT>(
      std::vector<size_t>({Nc,Nr}), level);
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
    T* outlowY = this->m_coeff->GetHalfTmpBuffPtr(0);
    T* outhighY = this->m_coeff->GetHalfTmpBuffPtr(1);
    for (int l=0; l<this->m_level; l++) {
      //Y filtering
      SeparableSubsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::f_l,
          typename WaveletSchemeT::f_h
          >::PerformSubsampledFilteringYRef(
        this->m_coeff->GetScaleShape(l).at(0),
        this->m_coeff->GetScaleShape(l).at(1),
        Accumulator<T,T,T,int>(inlow, outlowY,
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l).at(0)),
        Accumulator<T,T,T,int>(inlow, outhighY,
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l).at(0)));

      T* inlowY = outlowY;
      T* inhighY = outhighY;
      T* outlowYlowX;

      if (l+1==this->m_level) {
        outlowYlowX=this->m_coeff->GetLowSubspacePtr(l);
      } else {
        outlowYlowX=this->m_coeff->GetOutLowTmpBuffPtr();
      }

      //Now perform X filtering on lowpass Y
      SeparableSubsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::f_l,
          typename WaveletSchemeT::f_h
          >::PerformSubsampledFilteringXRef(
        this->m_coeff->GetScaleShape(l).at(0),
        this->m_coeff->GetScaleShape(l+1).at(0),
        this->m_coeff->GetScaleShape(l+1).at(1),
        Accumulator<T,T,T,int>(inlowY, outlowYlowX),
        Accumulator<T,T,T,int>(inlowY,
          this->m_coeff->GetHighSubspacePtr(l,0)));
      //Now perform X filtering on highpass Y
      SeparableSubsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::f_l,
          typename WaveletSchemeT::f_h
          >::PerformSubsampledFilteringXRef(
        this->m_coeff->GetScaleShape(l).at(0),
        this->m_coeff->GetScaleShape(l+1).at(0),
        this->m_coeff->GetScaleShape(l+1).at(1),
        Accumulator<T,T,T,int>(inhighY,
          this->m_coeff->GetHighSubspacePtr(l,1)),
        Accumulator<T,T,T,int>(inhighY,
          this->m_coeff->GetHighSubspacePtr(l,2)));

      //Update lowpass input and output, order is important here
      inlow=outlowYlowX;
    }
    return 1;
  }
  /// Backward wavelet transform: transpose of the forward transform
  virtual int backward() {

    for (int l=this->m_level; l>0; l--) {
      // Invert X lowpass/highpass filtering for lowpass Y
      SeparableUpsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::i_l,
          typename WaveletSchemeT::i_h
        >::PerformUpsampledFilteringXRef(
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l-1).at(0),
          this->m_coeff->GetScaleShape(l).at(1),
          this->m_coeff->GetHalfTmpBuffPtr(0),
          this->m_coeff->GetLowSubspacePtr(l-1),
          this->m_coeff->GetHighSubspacePtr(l-1,0));

      // Invert X lowpass/highpass filtering for highpass Y
      SeparableUpsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::i_l,
          typename WaveletSchemeT::i_h
        >::PerformUpsampledFilteringXRef(
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l-1).at(0),
          this->m_coeff->GetScaleShape(l).at(1),
          this->m_coeff->GetHalfTmpBuffPtr(1),
          this->m_coeff->GetHighSubspacePtr(l-1,1),
          this->m_coeff->GetHighSubspacePtr(l-1,2));

      //Update output buffer destination
      T* outlow;
      if (l<=1) {
        outlow=this->m_image;
      } else {
        outlow=this->m_coeff->GetLowSubspacePtr(l-2);
      }

      // Invert Y lowpass/highpass filtering
      SeparableUpsampledConvolutionEngine2D<T,
          typename WaveletSchemeT::i_l,
          typename WaveletSchemeT::i_h
        >::PerformUpsampledFilteringYRef(
          this->m_coeff->GetScaleShape(l-1).at(0),
          this->m_coeff->GetScaleShape(l).at(1),
          this->m_coeff->GetScaleShape(l-1).at(1),
          outlow,
          this->m_coeff->GetHalfTmpBuffPtr(0),
          this->m_coeff->GetHalfTmpBuffPtr(1));
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
using Dummy2_2D = Wavelet2D<T,PackedContainer2D<T>,Dummy2<T>>;
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

/** \class DTWavelet2D
 * \brief Inheritance of Wavelet class for the Dual Tree 2 dimensional case
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class CoeffContainerT, class DTWaveletSchemeT>
class DTWavelet2D : public Wavelet<T,CoeffContainerT, DTWaveletSchemeT> {
 public:
  /// Constructor with zero initialization
  DTWavelet2D()=default;

  /// Constructor : Wavelet from image
  DTWavelet2D(T* img, int Nc, int Nr, int Ns, bool doCycleSpinning,
      const std::string& wname, int level) : Wavelet<T,CoeffContainerT,
      DTWaveletSchemeT>(img, Nc, Nr, Ns, doCycleSpinning, wname, level) {
    this->m_coeff=std::make_unique<CoeffContainerT>(
      std::vector<size_t>({Nc,Nr}), level);
  }
  /// Default destructor
  virtual ~DTWavelet2D()=default;

  /// Forward wavelet tranform
  virtual int forward() {
    /**
     * At first step, input is simply input image
     * We chose to filter the most memory friendly direction at the end in
     * order to maximize performances as data size grows along filtering steps
     */
    T* inlow = this->m_image;
	T* outlowYReallowXReal;
	T* outlowYReallowXImag;
	T* outlowYImaglowXReal;
	T* outlowYImaglowXImag;
    int l = 0;

    //Building the quad tree
    if (this->m_level>=1) {
      // Divide the tree in two branch along the Y direction: Real and Imag
	  SeparableSubsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::f_l0r,
		  typename DTWaveletSchemeT::f_l0i,
		  typename DTWaveletSchemeT::f_h0r,
		  typename DTWaveletSchemeT::f_h0i
		  >::PerformSubsampledFilteringYRef(
		this->m_coeff->GetScaleShape(l).at(0),
		this->m_coeff->GetScaleShape(l).at(1),
		Accumulator<T,T,T,int>(inlow,
          this->m_coeff->GetHalfTmpBuffPtr(0,0),//LowYReal
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)),
		Accumulator<T,T,T,int>(inlow,
          this->m_coeff->GetHalfTmpBuffPtr(0,1),//LowYImag
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)),
		Accumulator<T,T,T,int>(inlow,
          this->m_coeff->GetHalfTmpBuffPtr(1,0),//HighYReal
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)),
		Accumulator<T,T,T,int>(inlow,
          this->m_coeff->GetHalfTmpBuffPtr(1,1),//HighYImag
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)));

	  if (l+1==this->m_level) {
		outlowYReallowXReal=this->m_coeff->GetLowSubspacePtr(l,0);
		outlowYReallowXImag=this->m_coeff->GetLowSubspacePtr(l,1);
		outlowYImaglowXReal=this->m_coeff->GetLowSubspacePtr(l,2);
		outlowYImaglowXImag=this->m_coeff->GetLowSubspacePtr(l,3);
	  } else {
		outlowYReallowXReal=this->m_coeff->GetOutLowTmpBuffPtr(0);
		outlowYReallowXImag=this->m_coeff->GetOutLowTmpBuffPtr(1);
		outlowYImaglowXReal=this->m_coeff->GetOutLowTmpBuffPtr(2);
		outlowYImaglowXImag=this->m_coeff->GetOutLowTmpBuffPtr(3);
	  }

	  //Now divide the tree in two branch along the X direction: Real in
      //First branch, first job: lowPassYReal -> YReal and YImag
	  SeparableSubsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::f_l0r,
		  typename DTWaveletSchemeT::f_l0i,
		  typename DTWaveletSchemeT::f_h0r,
		  typename DTWaveletSchemeT::f_h0i
		  >::PerformSubsampledFilteringXRef(
		this->m_coeff->GetScaleShape(l).at(0),
		this->m_coeff->GetScaleShape(l+1).at(0),
		this->m_coeff->GetScaleShape(l+1).at(1),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,0),//in: LowYReal
          outlowYReallowXReal),//out: LowYRealLowXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,0),//in: LowYReal
		  outlowYReallowXImag),//out: LowYRealLowXImag
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,0),//in: LowYReal
          this->m_coeff->GetHighSubspacePtr(l,0,0)),//out: LowYRealHighXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,0),//in: LowYReal
		  this->m_coeff->GetHighSubspacePtr(l,0,1)));//out: LowYRealHighXImag

      //First branch, second job: highPassYReal -> YReal and YImag
	  SeparableSubsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::f_l0r,
		  typename DTWaveletSchemeT::f_l0i,
		  typename DTWaveletSchemeT::f_h0r,
		  typename DTWaveletSchemeT::f_h0i
		  >::PerformSubsampledFilteringXRef(
		this->m_coeff->GetScaleShape(l).at(0),
		this->m_coeff->GetScaleShape(l+1).at(0),
		this->m_coeff->GetScaleShape(l+1).at(1),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,0),//in: HighYReal
          this->m_coeff->GetHighSubspacePtr(l,1,0)),//out: HighYRealLowXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,0),//in: HighYReal
		  this->m_coeff->GetHighSubspacePtr(l,1,1)),//out: HighYRealLowXImag
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,0),//in: HighYReal
          this->m_coeff->GetHighSubspacePtr(l,2,0)),//out: HighYRealHighXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,0),//in: HighYReal
		  this->m_coeff->GetHighSubspacePtr(l,2,1)));//out: HighYRealHighXImag

	  //Now divide the tree in two branch along the X direction: Imag in
      //Second branch, first job: lowPassYImag -> YReal and YImag
	  SeparableSubsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::f_l0r,
		  typename DTWaveletSchemeT::f_l0i,
		  typename DTWaveletSchemeT::f_h0r,
		  typename DTWaveletSchemeT::f_h0i
		  >::PerformSubsampledFilteringXRef(
		this->m_coeff->GetScaleShape(l).at(0),
		this->m_coeff->GetScaleShape(l+1).at(0),
		this->m_coeff->GetScaleShape(l+1).at(1),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,1),//in: LowYImag
          outlowYImaglowXReal),//out: LowYImagLowXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,1),//in: LowYImag
		  outlowYImaglowXImag),//out: LowYImagLowXImag
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,1),//in: LowYImag
          this->m_coeff->GetHighSubspacePtr(l,0,2)),//out: LowYRealHighXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,1),//in: LowYImag
		  this->m_coeff->GetHighSubspacePtr(l,0,3)));//out: LowYRealHighXImag

      //Second branch, second job: highPassYImag -> YReal and YImag
	  SeparableSubsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::f_l0r,
		  typename DTWaveletSchemeT::f_l0i,
		  typename DTWaveletSchemeT::f_h0r,
		  typename DTWaveletSchemeT::f_h0i
		  >::PerformSubsampledFilteringXRef(
		this->m_coeff->GetScaleShape(l).at(0),
		this->m_coeff->GetScaleShape(l+1).at(0),
		this->m_coeff->GetScaleShape(l+1).at(1),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,1),//in: HighYImag
          this->m_coeff->GetHighSubspacePtr(l,1,2)),//out: HighYRealLowXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,1),//in: HighYImag
		  this->m_coeff->GetHighSubspacePtr(l,1,3)),//out: HighYRealLowXImag
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,1),//in: HighYImag
          this->m_coeff->GetHighSubspacePtr(l,2,2)),//out: HighYRealHighXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,1),//in: HighYImag
		  this->m_coeff->GetHighSubspacePtr(l,2,3)));//out: HighYRealHighXImag
     }

    // Now each subtree can be processed separately
    for (int l=1; l<this->m_level; l++) {
	  if (l+1==this->m_level) {
		outlowYReallowXReal=this->m_coeff->GetLowSubspacePtr(l,0);
		outlowYReallowXImag=this->m_coeff->GetLowSubspacePtr(l,1);
		outlowYImaglowXReal=this->m_coeff->GetLowSubspacePtr(l,2);
		outlowYImaglowXImag=this->m_coeff->GetLowSubspacePtr(l,3);
	  } else {
		outlowYReallowXReal=this->m_coeff->GetOutLowTmpBuffPtr(0);
		outlowYReallowXImag=this->m_coeff->GetOutLowTmpBuffPtr(1);
		outlowYImaglowXReal=this->m_coeff->GetOutLowTmpBuffPtr(2);
		outlowYImaglowXImag=this->m_coeff->GetOutLowTmpBuffPtr(3);
	  }
      // Y filtering: Band 0 and 1: Real Y, Band 2 and 3: Imag Y
	  SeparableSubsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::f_lnr,
		  typename DTWaveletSchemeT::f_hnr,
		  typename DTWaveletSchemeT::f_lnr,
		  typename DTWaveletSchemeT::f_hnr,
		  typename DTWaveletSchemeT::f_lni,
		  typename DTWaveletSchemeT::f_hni,
		  typename DTWaveletSchemeT::f_lni,
		  typename DTWaveletSchemeT::f_hni
		  >::PerformSubsampledFilteringYRef(
		this->m_coeff->GetScaleShape(l).at(0),
		this->m_coeff->GetScaleShape(l).at(1),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetOutLowTmpBuffPtr(0),//in: lowYRealLowXReal
		  this->m_coeff->GetHalfTmpBuffPtr(0,0),//out: lowYReal
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetOutLowTmpBuffPtr(0),//in: lowYRealLowXReal
		  this->m_coeff->GetHalfTmpBuffPtr(1,0),//out: highYReal
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetOutLowTmpBuffPtr(1),//in: lowYRealLowXImag
		  this->m_coeff->GetHalfTmpBuffPtr(0,1),//out: lowYReal
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetOutLowTmpBuffPtr(1),//in: lowYRealLowXImag
		  this->m_coeff->GetHalfTmpBuffPtr(1,1),//out: highYReal
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetOutLowTmpBuffPtr(2),//in: lowYImagLowXReal
		  this->m_coeff->GetHalfTmpBuffPtr(0,2),//out: lowYImag
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetOutLowTmpBuffPtr(2),//in: lowYImagLowXReal
		  this->m_coeff->GetHalfTmpBuffPtr(1,2),//out: HighYImag
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetOutLowTmpBuffPtr(3),//in: lowYImagLowXImag
		  this->m_coeff->GetHalfTmpBuffPtr(0,3),//out: lowYImag
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetOutLowTmpBuffPtr(3),//in: lowYImagLowXImag
		  this->m_coeff->GetHalfTmpBuffPtr(1,3),//out: HighYImag
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(0)));

      //X filtering
	  SeparableSubsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::f_lnr,//Band0
		  typename DTWaveletSchemeT::f_hnr,
		  typename DTWaveletSchemeT::f_lnr,
		  typename DTWaveletSchemeT::f_hnr,
		  typename DTWaveletSchemeT::f_lni,//Band1
		  typename DTWaveletSchemeT::f_hni,
		  typename DTWaveletSchemeT::f_lni,
		  typename DTWaveletSchemeT::f_hni,
		  typename DTWaveletSchemeT::f_lnr,//Band2
		  typename DTWaveletSchemeT::f_hnr,
		  typename DTWaveletSchemeT::f_lnr,
		  typename DTWaveletSchemeT::f_hnr,
		  typename DTWaveletSchemeT::f_lni,//Band3
		  typename DTWaveletSchemeT::f_hni,
		  typename DTWaveletSchemeT::f_lni,
		  typename DTWaveletSchemeT::f_hni
		  >::PerformSubsampledFilteringXRef(
		this->m_coeff->GetScaleShape(l).at(0),
		this->m_coeff->GetScaleShape(l+1).at(0),
		this->m_coeff->GetScaleShape(l+1).at(1),
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,0),//in: LowYReal
          outlowYReallowXReal),//out: LowYRealLowXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,0),//in: LowYReal
		  this->m_coeff->GetHighSubspacePtr(l,0,0)),//out: LowYRealHighXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,0),//in: HighYReal
          this->m_coeff->GetHighSubspacePtr(l,1,0)),//out: HighYRealLowXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,0),//in: HighYReal
		  this->m_coeff->GetHighSubspacePtr(l,2,0)),//out: HighYRealHighXReal

		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,1),//in: LowYReal
          outlowYReallowXImag),//out: LowYRealLowXImag
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,1),//in: LowYReal
		  this->m_coeff->GetHighSubspacePtr(l,0,1)),//out: LowYRealHighXImag
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,1),//in: HighYReal
          this->m_coeff->GetHighSubspacePtr(l,1,1)),//out: HighYRealLowXImag
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,1),//in: HighYReal
		  this->m_coeff->GetHighSubspacePtr(l,2,1)),//out: HighYRealHighXImag

		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,2),//in: LowYImag
          outlowYImaglowXReal),//out: LowYImagLowXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,2),//in: LowYImag
		  this->m_coeff->GetHighSubspacePtr(l,0,2)),//out: LowYImagHighXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,2),//in: HighYImag
          this->m_coeff->GetHighSubspacePtr(l,1,2)),//out: HighYImagLowXReal
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,2),//in: HighYImag
		  this->m_coeff->GetHighSubspacePtr(l,2,2)),//out: HighYImagHighXReal

		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,3),//in: LowYImag
          outlowYImaglowXImag),//out: LowYImagLowXImag
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(0,3),//in: LowYImag
		  this->m_coeff->GetHighSubspacePtr(l,0,3)),//out: LowYImagHighXImag
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,3),//in: HighYImag
          this->m_coeff->GetHighSubspacePtr(l,1,3)),//out: HighYImagLowXImag
		Accumulator<T,T,T,int>(
          this->m_coeff->GetHalfTmpBuffPtr(1,3),//in: HighYImag
		  this->m_coeff->GetHighSubspacePtr(l,2,3)));//out: HighYImagHighXImag
    }


    // map the set of filtered signals to the real DTCWT mixture
    this->m_coeff->WaveletToCpx();

    return 1;
  }

  /// Backward wavelet transform: transpose of the forward transform
  virtual int backward() {

    // map the set of DTCWT coefficients to a simple set of filtered signals
    this->m_coeff->CpxToWavelet();

    for (int l=this->m_level; l>1; l--) {
      for (int bIdx : {0,2} ) {
		// Perform X filtering on band bIdx, low Y
		SeparableUpsampledConvolutionEngine2D<T,
			typename DTWaveletSchemeT::i_lnr,
			typename DTWaveletSchemeT::i_hnr
			>::PerformUpsampledFilteringXRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l-1).at(0),
		  this->m_coeff->GetScaleShape(l).at(1),
		  this->m_coeff->GetHalfTmpBuffPtr(0,bIdx),//out: LowYReal
		  this->m_coeff->GetLowSubspacePtr(l-1,bIdx),//in: LowRealYLowRealX
		  this->m_coeff->GetHighSubspacePtr(l-1,0,bIdx));//in:LowRealYHighRealX
		// Perform X filtering on band bIdx, High Y
		SeparableUpsampledConvolutionEngine2D<T,
			typename DTWaveletSchemeT::i_lnr,
			typename DTWaveletSchemeT::i_hnr
			>::PerformUpsampledFilteringXRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l-1).at(0),
		  this->m_coeff->GetScaleShape(l).at(1),
		  this->m_coeff->GetHalfTmpBuffPtr(1,bIdx),//out: HighYReal
		  this->m_coeff->GetHighSubspacePtr(l-1,1,bIdx),//in: HighRealYLowRealX
		  this->m_coeff->GetHighSubspacePtr(l-1,2,bIdx));//in:HighRealYHighRealX
      }
      for (int bIdx : {1,3} ) {
	    // Perform X filtering on band bIdx, low Y
		SeparableUpsampledConvolutionEngine2D<T,
			typename DTWaveletSchemeT::i_lni,
			typename DTWaveletSchemeT::i_hni
			>::PerformUpsampledFilteringXRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l-1).at(0),
		  this->m_coeff->GetScaleShape(l).at(1),
		  this->m_coeff->GetHalfTmpBuffPtr(0,bIdx),//out: LowYReal
		  this->m_coeff->GetLowSubspacePtr(l-1,bIdx),//in: LowRealYLowRealX
		  this->m_coeff->GetHighSubspacePtr(l-1,0,bIdx));//in:LowRealYHighRealX
		// Perform X filtering on band bIdx, High Y
		SeparableUpsampledConvolutionEngine2D<T,
			typename DTWaveletSchemeT::i_lni,
			typename DTWaveletSchemeT::i_hni
			>::PerformUpsampledFilteringXRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l-1).at(0),
		  this->m_coeff->GetScaleShape(l).at(1),
		  this->m_coeff->GetHalfTmpBuffPtr(1,bIdx),//out: HighYReal
		  this->m_coeff->GetHighSubspacePtr(l-1,1,bIdx),//in: HighRealYLowRealX
		  this->m_coeff->GetHighSubspacePtr(l-1,2,bIdx));//in:HighRealYHighRealX
      }
      //Now perform Y filtering
      for (int bIdx : {0,1} ) {
	    SeparableUpsampledConvolutionEngine2D<T,
		    typename DTWaveletSchemeT::i_lnr,
			typename DTWaveletSchemeT::i_hnr
			>::PerformUpsampledFilteringYRef(
		this->m_coeff->GetScaleShape(l-1).at(0),
		this->m_coeff->GetScaleShape(l).at(1),
		this->m_coeff->GetScaleShape(l-1).at(1),
		this->m_coeff->GetLowSubspacePtr(l-2,bIdx),//out: LowYLowX
		this->m_coeff->GetHalfTmpBuffPtr(0,bIdx), //in: LowY
		this->m_coeff->GetHalfTmpBuffPtr(1,bIdx));//in: HighY
      }
      for (int bIdx : {2,3} ) {
	    SeparableUpsampledConvolutionEngine2D<T,
		    typename DTWaveletSchemeT::i_lni,
			typename DTWaveletSchemeT::i_hni
			>::PerformUpsampledFilteringYRef(
		this->m_coeff->GetScaleShape(l-1).at(0),
		this->m_coeff->GetScaleShape(l).at(1),
		this->m_coeff->GetScaleShape(l-1).at(1),
		this->m_coeff->GetLowSubspacePtr(l-2,bIdx),//out: LowYLowX
		this->m_coeff->GetHalfTmpBuffPtr(0,bIdx), //in: LowY
		this->m_coeff->GetHalfTmpBuffPtr(1,bIdx));//in: HighY
      }
    }

    int l=1;
	T* lowYReallowXReal=this->m_coeff->GetLowSubspacePtr(l-1,0);
	T* lowYReallowXImag=this->m_coeff->GetLowSubspacePtr(l-1,1);
	T* lowYImaglowXReal=this->m_coeff->GetLowSubspacePtr(l-1,2);
	T* lowYImaglowXImag=this->m_coeff->GetLowSubspacePtr(l-1,3);
    
    if (this->m_level>=1) {

	  // Merge Band 0 and 1 in the X direction for lowYReal
	  SeparableUpsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::i_l0r,
		  typename DTWaveletSchemeT::i_l0i,
		  typename DTWaveletSchemeT::i_h0r,
		  typename DTWaveletSchemeT::i_h0i
		  >::PerformUpsampledFilteringXRef(
		this->m_coeff->GetScaleShape(l).at(0),
		this->m_coeff->GetScaleShape(l-1).at(0),
		this->m_coeff->GetScaleShape(l).at(1),
		this->m_coeff->GetHalfTmpBuffPtr(0,0),//LowYReal
		lowYReallowXReal,//in: LowRealYLowRealX
		lowYReallowXImag,//in: LowRealYLowImagX
		this->m_coeff->GetHighSubspacePtr(l-1,0,0),//in: LowRealYHighRealX
		this->m_coeff->GetHighSubspacePtr(l-1,0,1));//in: LowRealYHighImagX

      // Merge Band 0 and 1 in the X direction for highYReal
	  SeparableUpsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::i_l0r,
		  typename DTWaveletSchemeT::i_l0i,
		  typename DTWaveletSchemeT::i_h0r,
		  typename DTWaveletSchemeT::i_h0i
		  >::PerformUpsampledFilteringXRef(
		this->m_coeff->GetScaleShape(l).at(0),
		this->m_coeff->GetScaleShape(l-1).at(0),
		this->m_coeff->GetScaleShape(l).at(1),
		this->m_coeff->GetHalfTmpBuffPtr(1,0),//HighYReal
		this->m_coeff->GetHighSubspacePtr(l-1,1,0),//in: HighRealYLowRealX
		this->m_coeff->GetHighSubspacePtr(l-1,1,1),//in: HighRealYLowImagX
		this->m_coeff->GetHighSubspacePtr(l-1,2,0),//in: HighRealYHighRealX
		this->m_coeff->GetHighSubspacePtr(l-1,2,1));//in: HighRealYHighImagX

	  // Merge Band 2 and 3 in the X direction for lowYImag
	  SeparableUpsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::i_l0r,
		  typename DTWaveletSchemeT::i_l0i,
		  typename DTWaveletSchemeT::i_h0r,
		  typename DTWaveletSchemeT::i_h0i
		>::PerformUpsampledFilteringXRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l-1).at(0),
		  this->m_coeff->GetScaleShape(l).at(1),
          this->m_coeff->GetHalfTmpBuffPtr(0,1),//LowYImag
		  lowYImaglowXReal,//in: LowImagYLowRealX
		  lowYImaglowXImag,//in: LowImagYLowImagX
		  this->m_coeff->GetHighSubspacePtr(l-1,0,2),//in: LowImagYHighRealX
		  this->m_coeff->GetHighSubspacePtr(l-1,0,3));//in: LowImagYHighImagX

	  // Merge Band 2 and 3 in the X direction for highYImag
	  SeparableUpsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::i_l0r,
		  typename DTWaveletSchemeT::i_l0i,
		  typename DTWaveletSchemeT::i_h0r,
		  typename DTWaveletSchemeT::i_h0i
		>::PerformUpsampledFilteringXRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l-1).at(0),
          this->m_coeff->GetScaleShape(l).at(1),
	      this->m_coeff->GetHalfTmpBuffPtr(1,1),//HighYImag
		  this->m_coeff->GetHighSubspacePtr(l-1,1,2),//in: HighImagYLowRealX
		  this->m_coeff->GetHighSubspacePtr(l-1,1,3),//in: HighImagYLowImagX
		  this->m_coeff->GetHighSubspacePtr(l-1,2,2),//in: HighImagYHighRealX
		  this->m_coeff->GetHighSubspacePtr(l-1,2,3));//in: HighImagYHighImagX

      /****************************
       * Second Part: Y filtering
       ****************************/
      T* outlow = this->m_image;
	  
      // Now invert Y filtering, and merge everything
	  SeparableUpsampledConvolutionEngine2D<T,
		  typename DTWaveletSchemeT::i_l0r,
		  typename DTWaveletSchemeT::i_l0i,
		  typename DTWaveletSchemeT::i_h0r,
		  typename DTWaveletSchemeT::i_h0i
		>::PerformUpsampledFilteringYRef(
		  this->m_coeff->GetScaleShape(l-1).at(0),
		  this->m_coeff->GetScaleShape(l).at(1),
		  this->m_coeff->GetScaleShape(l-1).at(1),
		  outlow,
		  this->m_coeff->GetHalfTmpBuffPtr(0,0),//LowYReal
		  this->m_coeff->GetHalfTmpBuffPtr(0,1),//LowYImag
		  this->m_coeff->GetHalfTmpBuffPtr(1,0),//HighYReal
		  this->m_coeff->GetHalfTmpBuffPtr(1,1));//HighYImag
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
using PackedDTContainer2D =
  DTCoeffContainer2D<T,std::vector<T>>;

// Aliasing ugly types into more simple ones
template<typename T>
using dtwAnto97QSHIFT6_2D = 
  DTWavelet2D<T,PackedDTContainer2D<T>,dtwAnto97QSHIFT6<T>>;

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
 dtwAnto97QSHIFT6_2D<T> dtwAnto97QShift6_2D; 
};

#endif //WAVELET2D_H
