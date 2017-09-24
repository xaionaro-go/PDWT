#ifndef WAVELET3D_H
#define WAVELET3D_H

// Local
#include "Wavelet.h"

// STL
#include <vector>

// Local
#include "CoeffContainer.h"
#include "Filters.h"
#include "Separable.h"
#include "vectorization/Vectorization.h"


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
    this->m_coeff=std::make_unique<CoeffContainerT>(
      std::vector<size_t>({Nc,Nr,Ns}), level);
  }
  /// Default destructor
  virtual ~Wavelet3D()=default;

  /// Forward wavelet tranform
  virtual int forward() {
    /**
     * At first step, input is simply input image
     * We chose to filter the most memory friendly direction at the end in
     * order to mimic what may be done in overcomplete wavelet systems, where
     * data size grows along filtering steps
     */
    T* inlow = this->m_image;
    T* outlow = this->m_coeff->GetHalfTmpBuffPtr(2);
    // For every level
    for (int l=0; l<this->m_level; l++) {
      // Filter sequentially low and then high freq in z direction so that one
      // can use only one temporary buffer for the Z filtering stage
      for (int zFiltIdx=0; zFiltIdx<2; zFiltIdx++) {
        if (zFiltIdx==0) {
          //Z filtering low
          SeparableSubsampledConvolutionEngine3D<T,
            typename WaveletSchemeT::f_l
            >::PerformSubsampledFilteringZRef(
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l).at(1),
          this->m_coeff->GetScaleShape(l).at(2),
          Accumulator<T,T,T,int>(inlow,
            this->m_coeff->GetHalfTmpBuffPtr(0),
            this->m_coeff->GetScaleShape(l).at(0)*
              this->m_coeff->GetScaleShape(l).at(1),
            this->m_coeff->GetScaleShape(l).at(0)*
              this->m_coeff->GetScaleShape(l).at(1)));
        } else {
          //Z filtering high
          SeparableSubsampledConvolutionEngine3D<T,
            typename WaveletSchemeT::f_h
            >::PerformSubsampledFilteringZRef(
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l).at(1),
          this->m_coeff->GetScaleShape(l).at(2),
          Accumulator<T,T,T,int>(inlow,
            this->m_coeff->GetHalfTmpBuffPtr(0),
            this->m_coeff->GetScaleShape(l).at(0)*
              this->m_coeff->GetScaleShape(l).at(1),
            this->m_coeff->GetScaleShape(l).at(0)*
              this->m_coeff->GetScaleShape(l).at(1)));
            }

        for (int yFiltIdx=0; yFiltIdx<2; yFiltIdx++) {
          if (yFiltIdx==0) {
          //Y filtering low
          SeparableSubsampledConvolutionEngine3D<T,
            typename WaveletSchemeT::f_l
            >::PerformSubsampledFilteringYRef(
            this->m_coeff->GetScaleShape(l).at(0),
            this->m_coeff->GetScaleShape(l).at(1),
            this->m_coeff->GetScaleShape(l+1).at(1),
            this->m_coeff->GetScaleShape(l+1).at(2),
            Accumulator<T,T,T,int>(
              this->m_coeff->GetHalfTmpBuffPtr(0),
              this->m_coeff->GetHalfTmpBuffPtr(1),
              this->m_coeff->GetScaleShape(l).at(0),
              this->m_coeff->GetScaleShape(l).at(0)));
          } else {
            //Y filtering high
            SeparableSubsampledConvolutionEngine3D<T,
                typename WaveletSchemeT::f_h
                >::PerformSubsampledFilteringYRef(
              this->m_coeff->GetScaleShape(l).at(0),
              this->m_coeff->GetScaleShape(l).at(1),
              this->m_coeff->GetScaleShape(l+1).at(1),
              this->m_coeff->GetScaleShape(l+1).at(2),
              Accumulator<T,T,T,int>(
                this->m_coeff->GetHalfTmpBuffPtr(0),
                this->m_coeff->GetHalfTmpBuffPtr(1),
                this->m_coeff->GetScaleShape(l).at(0),
                this->m_coeff->GetScaleShape(l).at(0)));
          }

              auto sBandCalc = [=](auto xIdx){ return 
                this->m_coeff->GetHighSubspacePtr(l,zFiltIdx*4+yFiltIdx*2+xIdx-1);
              };
              T* outlowX; 
              // If this is the low freq projection
              if ((zFiltIdx==0)&&(yFiltIdx==0)) {
          if (l+1==this->m_level) {
            outlowX=this->m_coeff->GetLowSubspacePtr(l);
          } else {
            outlowX=outlow;
          }
              } else {
                outlowX=sBandCalc(0);
              }

          //Now perform X filtering
          SeparableSubsampledConvolutionEngine3D<T,
            typename WaveletSchemeT::f_l,
            typename WaveletSchemeT::f_h
            >::PerformSubsampledFilteringXRef(
          this->m_coeff->GetScaleShape(l).at(0),
          this->m_coeff->GetScaleShape(l+1).at(0),
          this->m_coeff->GetScaleShape(l+1).at(1),
          this->m_coeff->GetScaleShape(l+1).at(2),
          Accumulator<T,T,T,int>(this->m_coeff->GetHalfTmpBuffPtr(1),
                  outlowX),
          Accumulator<T,T,T,int>(this->m_coeff->GetHalfTmpBuffPtr(1),
            sBandCalc(1)));

        }
      }
      // Update address of low space projection
      if (l==0) {
        inlow = this->m_coeff->GetHalfTmpBuffPtr(2);
        outlow = this->m_coeff->GetHalfTmpBuffPtr(3);
      } else {
        std::swap(inlow,outlow);
      }
    }
    return 1;
  }


  /// Backward wavelet transform: transpose of the forward transpose
  virtual int backward() {
    int inTmpBuffIdx = (this->m_level%2==0) ? 3 : 2;
    int outTmpBuffIdx = (this->m_level%2==0) ? 2 : 3;
    T* inlow = this->m_coeff->GetHalfTmpBuffPtr(inTmpBuffIdx);
    T* outlow = this->m_coeff->GetHalfTmpBuffPtr(outTmpBuffIdx);
    for (int l=this->m_level; l>0; l--) {
      //Z low and high can be inverted independantly, and then added afterward
      for (int zFiltIdx=0; zFiltIdx<2; zFiltIdx++) {
        //Y low and high can be inverted sequentially, and then added
    for (int yFiltIdx=0; yFiltIdx<2; yFiltIdx++) {
          // returns ptr to coefficients for the current (x,y,z) configuration
          auto sBandCalc = [=](auto xIdx){ return 
            this->m_coeff->GetHighSubspacePtr(l-1,zFiltIdx*4+yFiltIdx*2+xIdx-1);
          };
          T* inlowX;
          // If this is the low freq projection
          if ((zFiltIdx==0)&&(yFiltIdx==0)) {
            if (l==this->m_level) {
              inlowX=this->m_coeff->GetLowSubspacePtr(l-1);
            } else {
              inlowX=inlow;
            }
          } else {
            inlowX=sBandCalc(0);
          }

          // Invert X lowpass/highpass filtering for any Y/Z combination
        SeparableUpsampledConvolutionEngine3D<T,
        SubsampledAccumulator,
        typename WaveletSchemeT::i_l,
        typename WaveletSchemeT::i_h
        >::PerformUpsampledFilteringXRef(
          this->m_coeff->GetScaleShape(l).at(0),
      this->m_coeff->GetScaleShape(l-1).at(0),
      this->m_coeff->GetScaleShape(l).at(1),
      this->m_coeff->GetScaleShape(l).at(2),
      this->m_coeff->GetHalfTmpBuffPtr(1),
      inlowX,
      sBandCalc(1));

          if (yFiltIdx==0) {
      // Invert Y lowpass filtering only
      SeparableUpsampledConvolutionEngine3D<T,
                SubsampledAccumulator,
        typename WaveletSchemeT::i_l
        >::PerformUpsampledFilteringYRef(
        this->m_coeff->GetScaleShape(l-1).at(0),
        this->m_coeff->GetScaleShape(l).at(1),
        this->m_coeff->GetScaleShape(l-1).at(1),
        this->m_coeff->GetScaleShape(l).at(2),
        this->m_coeff->GetHalfTmpBuffPtr(0),
        this->m_coeff->GetHalfTmpBuffPtr(1)); 
          } else { //Perform update instead of write, see Accumulator type
      // Invert Y highpass filtering only
      SeparableUpsampledConvolutionEngine3D<T,
                SubsampledAccumulatorUpdate,
        typename WaveletSchemeT::i_h
        >::PerformUpsampledFilteringYRef(
        this->m_coeff->GetScaleShape(l-1).at(0),
        this->m_coeff->GetScaleShape(l).at(1),
        this->m_coeff->GetScaleShape(l-1).at(1),
        this->m_coeff->GetScaleShape(l).at(2),
        this->m_coeff->GetHalfTmpBuffPtr(0),
        this->m_coeff->GetHalfTmpBuffPtr(1));
          }
        }

    T* outlowZ;
    if (l<=1) {
      outlowZ=this->m_image;
    } else {
      outlowZ=outlow;
    }

        if (zFiltIdx==0) {
      // Invert Z lowpass filtering
      SeparableUpsampledConvolutionEngine3D<T,
              SubsampledAccumulator,
        typename WaveletSchemeT::i_l
      >::PerformUpsampledFilteringZRef(
        this->m_coeff->GetScaleShape(l-1).at(0),
        this->m_coeff->GetScaleShape(l-1).at(1),
        this->m_coeff->GetScaleShape(l).at(2),
        this->m_coeff->GetScaleShape(l-1).at(2),
        outlowZ,
        this->m_coeff->GetHalfTmpBuffPtr(0));
        } else { //Perform update instead of write, see Accumulator tyoe
      // Invert Z highpass filtering
      SeparableUpsampledConvolutionEngine3D<T,
              SubsampledAccumulatorUpdate,
        typename WaveletSchemeT::i_h
      >::PerformUpsampledFilteringZRef(
        this->m_coeff->GetScaleShape(l-1).at(0),
        this->m_coeff->GetScaleShape(l-1).at(1),
        this->m_coeff->GetScaleShape(l).at(2),
        this->m_coeff->GetScaleShape(l-1).at(2),
        outlowZ,
        this->m_coeff->GetHalfTmpBuffPtr(0));
        }
    }
      std::swap(inlow, outlow);
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
using PackedContainer3D =
  CoeffContainer3D<T,std::vector<T,PackAllocator<T>>>;

// Aliasing ugly types into more simple ones
template<typename T>
using Daub2_3D = Wavelet3D<T,PackedContainer3D<T>,Daub2<T>>;
template<typename T>
using Daub3_3D = Wavelet3D<T,PackedContainer3D<T>,Daub3<T>>;
template<typename T>
using Daub4_3D = Wavelet3D<T,PackedContainer3D<T>,Daub4<T>>;
template<typename T>
using Daub5_3D = Wavelet3D<T,PackedContainer3D<T>,Daub5<T>>;
template<typename T>
using Anto97_BiOrth_3D = Wavelet3D<T,PackedContainer3D<T>,Anto97_BiOrth<T>>;
template<typename T>
using QSHIFT6_Orth_3D = Wavelet3D<T,PackedContainer3D<T>,QSHIFT6_Orth<T>>;
template<typename T>
using REVERSE_QSHIFT6_Orth_3D = 
  Wavelet3D<T,PackedContainer3D<T>,REVERSE_QSHIFT6_Orth<T>>;

/** \class DTWavelet3D
 * \brief Inheritance of Wavelet class for the Dual Tree 2 dimensional case
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class CoeffContainerT, class DTWaveletSchemeT>
class DTWavelet3D : public Wavelet<T,CoeffContainerT, DTWaveletSchemeT> {
 public:
  /// Constructor with zero initialization
  DTWavelet3D()=default;

  /// Constructor : Wavelet from image
  DTWavelet3D(T* img, int Nc, int Nr, int Ns, bool doCycleSpinning,
      const std::string& wname, int level) : Wavelet<T,CoeffContainerT,
      DTWaveletSchemeT>(img, Nc, Nr, Ns, doCycleSpinning, wname, level) {
    this->m_coeff=std::make_unique<CoeffContainerT>(
      std::vector<size_t>({Nc,Nr,Ns}), level);
  }
  /// Default destructor
  virtual ~DTWavelet3D()=default;

  /// Forward wavelet tranform
  virtual int forward() {
    /**
     * At first step, input is simply input image
     * We chose to filter the most memory friendly direction at the end in
     * order to maximize performances as data size grows along filtering steps
     */
    T* inlow = this->m_image;
    int l = 0;

    //First build the octree, then perform regular CWT on each tree
    if (this->m_level>=1) {
      /********************
       * Z Filtering step *
       ********************/
      //Filter with first zReal filter, and then zImag filter
      for (int zImagStatus : {0,1}) {
        // Filter sequentially low and then high freq in z direction so that one
        // can use only one temporary buffer for the Z filtering stage
        for (int zFiltIdx=0; zFiltIdx<2; zFiltIdx++) {
          if (zImagStatus==0) {
            //Z Real Filters
            if (zFiltIdx==0) {
              //Z Real Filters low
              SeparableSubsampledConvolutionEngine3D<T,
                  typename DTWaveletSchemeT::f_l0r
                  >::PerformSubsampledFilteringZRef(
                this->m_coeff->GetScaleShape(l).at(0),
                this->m_coeff->GetScaleShape(l).at(1),
                this->m_coeff->GetScaleShape(l).at(2),
                Accumulator<T,T,T,int>(inlow,
                  this->m_coeff->GetHalfTmpBuffPtr(0),
                  this->m_coeff->GetScaleShape(l).at(0)*
                    this->m_coeff->GetScaleShape(l).at(1),
                  this->m_coeff->GetScaleShape(l).at(0)*
                    this->m_coeff->GetScaleShape(l).at(1)));
            } else {
              //Z Real Filters high
              SeparableSubsampledConvolutionEngine3D<T,
                  typename DTWaveletSchemeT::f_h0r
                  >::PerformSubsampledFilteringZRef(
                this->m_coeff->GetScaleShape(l).at(0),
                this->m_coeff->GetScaleShape(l).at(1),
                this->m_coeff->GetScaleShape(l).at(2),
                Accumulator<T,T,T,int>(inlow,
                  this->m_coeff->GetHalfTmpBuffPtr(0),
                  this->m_coeff->GetScaleShape(l).at(0)*
                  this->m_coeff->GetScaleShape(l).at(1),
                  this->m_coeff->GetScaleShape(l).at(0)*
                  this->m_coeff->GetScaleShape(l).at(1)));
        }
          } else {
            //Z Imag Filters
          if (zFiltIdx==0) {
              //Z Imag Filters low
              SeparableSubsampledConvolutionEngine3D<T,
                  typename DTWaveletSchemeT::f_l0i
                  >::PerformSubsampledFilteringZRef(
                this->m_coeff->GetScaleShape(l).at(0),
                this->m_coeff->GetScaleShape(l).at(1),
                this->m_coeff->GetScaleShape(l).at(2),
                Accumulator<T,T,T,int>(inlow,
                  this->m_coeff->GetHalfTmpBuffPtr(0),
                  this->m_coeff->GetScaleShape(l).at(0)*
                  this->m_coeff->GetScaleShape(l).at(1),
                  this->m_coeff->GetScaleShape(l).at(0)*
                  this->m_coeff->GetScaleShape(l).at(1)));
        } else {
              //Z Imag Filters high
              SeparableSubsampledConvolutionEngine3D<T,
                  typename DTWaveletSchemeT::f_h0i
                  >::PerformSubsampledFilteringZRef(
                this->m_coeff->GetScaleShape(l).at(0),
                this->m_coeff->GetScaleShape(l).at(1),
                this->m_coeff->GetScaleShape(l).at(2),
                Accumulator<T,T,T,int>(inlow,
                  this->m_coeff->GetHalfTmpBuffPtr(0),
                  this->m_coeff->GetScaleShape(l).at(0)*
                  this->m_coeff->GetScaleShape(l).at(1),
                  this->m_coeff->GetScaleShape(l).at(0)*
                  this->m_coeff->GetScaleShape(l).at(1)));
        }
          }
          /********************
           * Y Filtering step *
           ********************/
	  //Filter with first yReal filter, and then yImag filter
	  for (int yImagStatus : {0,1}) {
	    // Filter sequentially low/high freq in z direction so that one
	    // can use only one temporary buffer for the Z filtering stage
	    for (int yFiltIdx=0; yFiltIdx<2; yFiltIdx++) {
	      if (yImagStatus==0) {
		//Y Real Filters
		if (yFiltIdx==0) {
		  //Y Real Filters low
		  SeparableSubsampledConvolutionEngine3D<T,
		      typename DTWaveletSchemeT::f_l0r
		      >::PerformSubsampledFilteringYRef(
		    this->m_coeff->GetScaleShape(l).at(0),
		    this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l+1).at(1),
		    this->m_coeff->GetScaleShape(l+1).at(2),
		    Accumulator<T,T,T,int>(
		      this->m_coeff->GetHalfTmpBuffPtr(0),
		      this->m_coeff->GetHalfTmpBuffPtr(1),
		      this->m_coeff->GetScaleShape(l).at(0),
		      this->m_coeff->GetScaleShape(l).at(0)));
		} else { 
		  //Y Real Filters high
		  SeparableSubsampledConvolutionEngine3D<T,
		      typename DTWaveletSchemeT::f_h0r
		      >::PerformSubsampledFilteringYRef(
		    this->m_coeff->GetScaleShape(l).at(0),
		    this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l+1).at(1),
		    this->m_coeff->GetScaleShape(l+1).at(2),
		    Accumulator<T,T,T,int>(
		      this->m_coeff->GetHalfTmpBuffPtr(0),
		      this->m_coeff->GetHalfTmpBuffPtr(1),
		      this->m_coeff->GetScaleShape(l).at(0),
		      this->m_coeff->GetScaleShape(l).at(0)));
		}
	      } else {
		//Y Imag Filters
		if (yFiltIdx==0) {
		  //Y Imag Filters low
		  SeparableSubsampledConvolutionEngine3D<T,
		      typename DTWaveletSchemeT::f_l0i
		      >::PerformSubsampledFilteringYRef(
		    this->m_coeff->GetScaleShape(l).at(0),
		    this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l+1).at(1),
		    this->m_coeff->GetScaleShape(l+1).at(2),
		    Accumulator<T,T,T,int>(
		      this->m_coeff->GetHalfTmpBuffPtr(0),
		      this->m_coeff->GetHalfTmpBuffPtr(1),
		      this->m_coeff->GetScaleShape(l).at(0),
		      this->m_coeff->GetScaleShape(l).at(0)));
		} else {
		  //Y Imag Filters high
		  SeparableSubsampledConvolutionEngine3D<T,
		      typename DTWaveletSchemeT::f_h0i
		      >::PerformSubsampledFilteringYRef(
		    this->m_coeff->GetScaleShape(l).at(0),
		    this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l+1).at(1),
		    this->m_coeff->GetScaleShape(l+1).at(2),
		    Accumulator<T,T,T,int>(
		      this->m_coeff->GetHalfTmpBuffPtr(0),
		      this->m_coeff->GetHalfTmpBuffPtr(1),
		      this->m_coeff->GetScaleShape(l).at(0),
		      this->m_coeff->GetScaleShape(l).at(0)));
		}
	      }
              /********************
               * X Filtering step *
               ********************/
              T *outlowXReal, *outlowXImag;
              int bandIdxReal=4*zImagStatus+2*yImagStatus;
              int bandIdxImag=bandIdxReal+1;
              auto sBandCalc = [=](auto xIdx, auto bandIdx){
                return this->m_coeff->GetHighSubspacePtr(
                  l ,zFiltIdx*4+yFiltIdx*2+xIdx-1, bandIdx);
              };
              if ((zFiltIdx==0)&&(yFiltIdx==0)) {
                if (l+1==this->m_level) {
                  outlowXReal=this->m_coeff->GetLowSubspacePtr(l,bandIdxReal);
                  outlowXImag=this->m_coeff->GetLowSubspacePtr(l,bandIdxImag);
                } else {
                  outlowXReal=this->m_coeff->GetOutLowTmpBuffPtr(bandIdxReal);
                  outlowXImag=this->m_coeff->GetOutLowTmpBuffPtr(bandIdxImag);
                }
              } else {
                outlowXReal=sBandCalc(0, bandIdxReal);
                outlowXImag=sBandCalc(0, bandIdxImag);
              }
              //Perform full filtering in once
              SeparableSubsampledConvolutionEngine3D<T,
                  typename DTWaveletSchemeT::f_l0r,
                  typename DTWaveletSchemeT::f_h0r,
                  typename DTWaveletSchemeT::f_l0i,
                  typename DTWaveletSchemeT::f_h0i
                  >::PerformSubsampledFilteringXRef(
                this->m_coeff->GetScaleShape(l).at(0),
                this->m_coeff->GetScaleShape(l+1).at(0),
                this->m_coeff->GetScaleShape(l+1).at(1),
                this->m_coeff->GetScaleShape(l+1).at(2),
                Accumulator<T,T,T,int>(this->m_coeff->GetHalfTmpBuffPtr(1),
                  outlowXReal),
                Accumulator<T,T,T,int>(this->m_coeff->GetHalfTmpBuffPtr(1),
                  sBandCalc(1, bandIdxReal)),
                Accumulator<T,T,T,int>(this->m_coeff->GetHalfTmpBuffPtr(1),
                  outlowXImag),
                Accumulator<T,T,T,int>(this->m_coeff->GetHalfTmpBuffPtr(1),
                  sBandCalc(1, bandIdxImag)));
            }
          }
        }

        // We process each subtree sequentially
        for (int bandIdx=0; bandIdx<8; bandIdx++) {
          T* inlow = this->m_coeff->GetOutLowTmpBuffPtr(bandIdx);
          T* outlow = this->m_coeff->GetHalfTmpBuffPtr(3, bandIdx);
          for (int l=1; l<this->m_level; l++) {
	    // Filter sequentially low and then high freq in z direction so that one
	    // can use only one temporary buffer for the Z filtering stage
	    for (int zFiltIdx=0; zFiltIdx<2; zFiltIdx++) {
	      if (bandIdx<4) {
		if (zFiltIdx==0) {
		  //Z filtering low
		  SeparableSubsampledConvolutionEngine3D<T,
		    typename DTWaveletSchemeT::f_lnr
		    >::PerformSubsampledFilteringZRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(1),
		  this->m_coeff->GetScaleShape(l).at(2),
		  Accumulator<T,T,T,int>(inlow,
		    this->m_coeff->GetHalfTmpBuffPtr(0),
		    this->m_coeff->GetScaleShape(l).at(0)*
		      this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l).at(0)*
		      this->m_coeff->GetScaleShape(l).at(1)));
		} else {
		  //Z filtering high
		  SeparableSubsampledConvolutionEngine3D<T,
		    typename DTWaveletSchemeT::f_hnr
		    >::PerformSubsampledFilteringZRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(1),
		  this->m_coeff->GetScaleShape(l).at(2),
		  Accumulator<T,T,T,int>(inlow,
		    this->m_coeff->GetHalfTmpBuffPtr(0),
		    this->m_coeff->GetScaleShape(l).at(0)*
		      this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l).at(0)*
		      this->m_coeff->GetScaleShape(l).at(1)));
		 }
	       } else {
		if (zFiltIdx==0) {
		  //Z filtering low
		  SeparableSubsampledConvolutionEngine3D<T,
		    typename DTWaveletSchemeT::f_lni
		    >::PerformSubsampledFilteringZRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(1),
		  this->m_coeff->GetScaleShape(l).at(2),
		  Accumulator<T,T,T,int>(inlow,
		    this->m_coeff->GetHalfTmpBuffPtr(0),
		    this->m_coeff->GetScaleShape(l).at(0)*
		      this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l).at(0)*
		      this->m_coeff->GetScaleShape(l).at(1)));
		} else {
		  //Z filtering high
		  SeparableSubsampledConvolutionEngine3D<T,
		    typename DTWaveletSchemeT::f_hni
		    >::PerformSubsampledFilteringZRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l).at(1),
		  this->m_coeff->GetScaleShape(l).at(2),
		  Accumulator<T,T,T,int>(inlow,
		    this->m_coeff->GetHalfTmpBuffPtr(0),
		    this->m_coeff->GetScaleShape(l).at(0)*
		      this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l).at(0)*
		      this->m_coeff->GetScaleShape(l).at(1)));
	        }
              }

	      for (int yFiltIdx=0; yFiltIdx<2; yFiltIdx++) {
		if ((bandIdx%4)<2) {
		  if (yFiltIdx==0) {
		    //Y filtering low
		    SeparableSubsampledConvolutionEngine3D<T,
		      typename DTWaveletSchemeT::f_lnr
		      >::PerformSubsampledFilteringYRef(
		      this->m_coeff->GetScaleShape(l).at(0),
		      this->m_coeff->GetScaleShape(l).at(1),
		      this->m_coeff->GetScaleShape(l+1).at(1),
		      this->m_coeff->GetScaleShape(l+1).at(2),
		      Accumulator<T,T,T,int>(
			this->m_coeff->GetHalfTmpBuffPtr(0),
			this->m_coeff->GetHalfTmpBuffPtr(1),
			this->m_coeff->GetScaleShape(l).at(0),
			this->m_coeff->GetScaleShape(l).at(0)));
		  } else {
		    //Y filtering high
		    SeparableSubsampledConvolutionEngine3D<T,
			typename DTWaveletSchemeT::f_hnr
			>::PerformSubsampledFilteringYRef(
		      this->m_coeff->GetScaleShape(l).at(0),
		      this->m_coeff->GetScaleShape(l).at(1),
		      this->m_coeff->GetScaleShape(l+1).at(1),
		      this->m_coeff->GetScaleShape(l+1).at(2),
		      Accumulator<T,T,T,int>(
			this->m_coeff->GetHalfTmpBuffPtr(0),
			this->m_coeff->GetHalfTmpBuffPtr(1),
			this->m_coeff->GetScaleShape(l).at(0),
			this->m_coeff->GetScaleShape(l).at(0)));
		  }
		} else {
		  if (yFiltIdx==0) {
		    //Y filtering low
		    SeparableSubsampledConvolutionEngine3D<T,
		      typename DTWaveletSchemeT::f_lni
		      >::PerformSubsampledFilteringYRef(
		      this->m_coeff->GetScaleShape(l).at(0),
		      this->m_coeff->GetScaleShape(l).at(1),
		      this->m_coeff->GetScaleShape(l+1).at(1),
		      this->m_coeff->GetScaleShape(l+1).at(2),
		      Accumulator<T,T,T,int>(
			this->m_coeff->GetHalfTmpBuffPtr(0),
			this->m_coeff->GetHalfTmpBuffPtr(1),
			this->m_coeff->GetScaleShape(l).at(0),
			this->m_coeff->GetScaleShape(l).at(0)));
		  } else {
		    //Y filtering high
		    SeparableSubsampledConvolutionEngine3D<T,
			typename DTWaveletSchemeT::f_hni
			>::PerformSubsampledFilteringYRef(
		      this->m_coeff->GetScaleShape(l).at(0),
		      this->m_coeff->GetScaleShape(l).at(1),
		      this->m_coeff->GetScaleShape(l+1).at(1),
		      this->m_coeff->GetScaleShape(l+1).at(2),
		      Accumulator<T,T,T,int>(
			this->m_coeff->GetHalfTmpBuffPtr(0),
			this->m_coeff->GetHalfTmpBuffPtr(1),
			this->m_coeff->GetScaleShape(l).at(0),
			this->m_coeff->GetScaleShape(l).at(0)));
		  }

		}

		auto sBandCalc = [=](auto xIdx, auto band){ return 
		  this->m_coeff->GetHighSubspacePtr(
		    l,zFiltIdx*4+yFiltIdx*2+xIdx-1, band);
		};
		T* outlowX;
		// If this is the low freq projection
		if ((zFiltIdx==0)&&(yFiltIdx==0)) {
		  if (l+1==this->m_level) {
		    outlowX=this->m_coeff->GetLowSubspacePtr(l, bandIdx);
		  } else {
		    outlowX=outlow;
		  }
		} else {
		  outlowX=sBandCalc(0, bandIdx);
		}

		if ((bandIdx%2)==0) {
		  //Now perform X filtering
		  SeparableSubsampledConvolutionEngine3D<T,
		    typename DTWaveletSchemeT::f_lnr,
		    typename DTWaveletSchemeT::f_hnr
		    >::PerformSubsampledFilteringXRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l+1).at(0),
		  this->m_coeff->GetScaleShape(l+1).at(1),
		  this->m_coeff->GetScaleShape(l+1).at(2),
		  Accumulator<T,T,T,int>(this->m_coeff->GetHalfTmpBuffPtr(1),
		    outlowX),
		  Accumulator<T,T,T,int>(this->m_coeff->GetHalfTmpBuffPtr(1),
		    sBandCalc(1, bandIdx)));
		} else {
		  //Now perform X filtering
		  SeparableSubsampledConvolutionEngine3D<T,
		    typename DTWaveletSchemeT::f_lni,
		    typename DTWaveletSchemeT::f_hni
		    >::PerformSubsampledFilteringXRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l+1).at(0),
		  this->m_coeff->GetScaleShape(l+1).at(1),
		  this->m_coeff->GetScaleShape(l+1).at(2),
		  Accumulator<T,T,T,int>(this->m_coeff->GetHalfTmpBuffPtr(1),
		    outlowX),
		  Accumulator<T,T,T,int>(this->m_coeff->GetHalfTmpBuffPtr(1),
		    sBandCalc(1, bandIdx)));

		}
	      }
	    }
	    // Update address of low space projection
	    std::swap(inlow,outlow);
          }
        }
      }

      // map the set of filtered signals to the real DTCWT mixture
      this->m_coeff->WaveletToCpx();
    }
    return 1;
  }

  /// Backward wavelet transform: transpose of the forward transform
  virtual int backward() {

    // map the set of DTCWT coefficients to a simple set of filtered signals
    this->m_coeff->CpxToWavelet();

    if (this->m_level>=1) {
      //Every band can be reconstructed independantly
      for (int bandIdx=0; bandIdx<8; bandIdx++) {
        int inTmpBuffIdx = (this->m_level%2==0) ? 3 : 2;
        int outTmpBuffIdx = (this->m_level%2==0) ? 2 : 3;
        T* inlow = this->m_coeff->GetHalfTmpBuffPtr(inTmpBuffIdx, bandIdx);
        T* outlow = this->m_coeff->GetHalfTmpBuffPtr(outTmpBuffIdx, bandIdx);
        for (int l=this->m_level; l>1; l--) {
          //Z low and high can be inverted independantly, and then added after
          for (int zFiltIdx=0; zFiltIdx<2; zFiltIdx++) {
            //Y low and high can be inverted sequentially, and then added
            for (int yFiltIdx=0; yFiltIdx<2; yFiltIdx++) {
              // returns ptr to coefficients for the current (x,y,z) confg
              auto sBandCalc = [=](auto xIdx){ return 
                this->m_coeff->GetHighSubspacePtr(
                  l-1,zFiltIdx*4+yFiltIdx*2+xIdx-1, bandIdx);
              };
              T* inlowX;
              // If this is the low freq projection
              if ((zFiltIdx==0)&&(yFiltIdx==0)) {
                if (l==this->m_level) {
                  inlowX=this->m_coeff->GetLowSubspacePtr(
                    this->m_level,bandIdx);
                } else {
                  inlowX=inlow;
                }
              } else {
                inlowX=sBandCalc(0);
              }

              // Invert X lowpass/highpass filtering for any Y/Z combination
              if (bandIdx%2==0) {
                SeparableUpsampledConvolutionEngine3D<T,
                    SubsampledAccumulator,
                    typename DTWaveletSchemeT::i_lnr,
                    typename DTWaveletSchemeT::i_hnr
                  >::PerformUpsampledFilteringXRef(
                    this->m_coeff->GetScaleShape(l).at(0),
                  this->m_coeff->GetScaleShape(l-1).at(0),
                  this->m_coeff->GetScaleShape(l).at(1),
                  this->m_coeff->GetScaleShape(l).at(2),
                  this->m_coeff->GetHalfTmpBuffPtr(1),
                  inlowX,
                  sBandCalc(1));
              } else {
                SeparableUpsampledConvolutionEngine3D<T,
                    SubsampledAccumulator,
                    typename DTWaveletSchemeT::i_lni,
                    typename DTWaveletSchemeT::i_hni
                  >::PerformUpsampledFilteringXRef(
                    this->m_coeff->GetScaleShape(l).at(0),
                  this->m_coeff->GetScaleShape(l-1).at(0),
                  this->m_coeff->GetScaleShape(l).at(1),
                  this->m_coeff->GetScaleShape(l).at(2),
                  this->m_coeff->GetHalfTmpBuffPtr(1),
                  inlowX,
                  sBandCalc(1));
              }
              // Then filter along Y axis
              if ((bandIdx%4)<2) {
                if (yFiltIdx==0) {
                  // Invert Y lowpass filtering only
                  SeparableUpsampledConvolutionEngine3D<T,
                      SubsampledAccumulator,
                      typename DTWaveletSchemeT::i_lnr
                    >::PerformUpsampledFilteringYRef(
                    this->m_coeff->GetScaleShape(l-1).at(0),
                    this->m_coeff->GetScaleShape(l).at(1),
                    this->m_coeff->GetScaleShape(l-1).at(1),
                    this->m_coeff->GetScaleShape(l).at(2),
                    this->m_coeff->GetHalfTmpBuffPtr(0),
                    this->m_coeff->GetHalfTmpBuffPtr(1)); 
                } else {
                  // Invert Y highpass filtering only
                  SeparableUpsampledConvolutionEngine3D<T,
                      SubsampledAccumulatorUpdate,
                      typename DTWaveletSchemeT::i_hnr
                    >::PerformUpsampledFilteringYRef(
                    this->m_coeff->GetScaleShape(l-1).at(0),
                    this->m_coeff->GetScaleShape(l).at(1),
                    this->m_coeff->GetScaleShape(l-1).at(1),
                    this->m_coeff->GetScaleShape(l).at(2),
                    this->m_coeff->GetHalfTmpBuffPtr(0),
                    this->m_coeff->GetHalfTmpBuffPtr(1));
                }
              } else {
                if (yFiltIdx==0) {
                  // Invert Y lowpass filtering only
                  SeparableUpsampledConvolutionEngine3D<T,
                      SubsampledAccumulator,
                      typename DTWaveletSchemeT::i_lni
                    >::PerformUpsampledFilteringYRef(
                    this->m_coeff->GetScaleShape(l-1).at(0),
                    this->m_coeff->GetScaleShape(l).at(1),
                    this->m_coeff->GetScaleShape(l-1).at(1),
                    this->m_coeff->GetScaleShape(l).at(2),
                    this->m_coeff->GetHalfTmpBuffPtr(0),
                    this->m_coeff->GetHalfTmpBuffPtr(1)); 
                } else {
                  // Invert Y highpass filtering only
                  SeparableUpsampledConvolutionEngine3D<T,
                      SubsampledAccumulatorUpdate,
                      typename DTWaveletSchemeT::i_hni
                    >::PerformUpsampledFilteringYRef(
                    this->m_coeff->GetScaleShape(l-1).at(0),
                    this->m_coeff->GetScaleShape(l).at(1),
                    this->m_coeff->GetScaleShape(l-1).at(1),
                    this->m_coeff->GetScaleShape(l).at(2),
                    this->m_coeff->GetHalfTmpBuffPtr(0),
                    this->m_coeff->GetHalfTmpBuffPtr(1));
                }
              }
            }
  
            if (bandIdx<4) {
              if (zFiltIdx==0) {
                // Invert Z lowpass filtering
                SeparableUpsampledConvolutionEngine3D<T,
                    SubsampledAccumulator,
                    typename DTWaveletSchemeT::i_lnr
                  >::PerformUpsampledFilteringZRef(
                  this->m_coeff->GetScaleShape(l-1).at(0),
                  this->m_coeff->GetScaleShape(l-1).at(1),
                  this->m_coeff->GetScaleShape(l).at(2),
                  this->m_coeff->GetScaleShape(l-1).at(2),
                  outlow,
                  this->m_coeff->GetHalfTmpBuffPtr(0));
              } else { //Perform update instead of write, see Accumulator tyoe
                // Invert Z highpass filtering
                SeparableUpsampledConvolutionEngine3D<T,
                    SubsampledAccumulatorUpdate,
                    typename DTWaveletSchemeT::i_hnr
                  >::PerformUpsampledFilteringZRef(
                  this->m_coeff->GetScaleShape(l-1).at(0),
                  this->m_coeff->GetScaleShape(l-1).at(1),
                  this->m_coeff->GetScaleShape(l).at(2),
                  this->m_coeff->GetScaleShape(l-1).at(2),
                  outlow,
                  this->m_coeff->GetHalfTmpBuffPtr(0));
              }
            } else {
              if (zFiltIdx==0) {
                // Invert Z lowpass filtering
                SeparableUpsampledConvolutionEngine3D<T,
                    SubsampledAccumulator,
                    typename DTWaveletSchemeT::i_lni
                  >::PerformUpsampledFilteringZRef(
                  this->m_coeff->GetScaleShape(l-1).at(0),
                  this->m_coeff->GetScaleShape(l-1).at(1),
                  this->m_coeff->GetScaleShape(l).at(2),
                  this->m_coeff->GetScaleShape(l-1).at(2),
                  outlow,
                  this->m_coeff->GetHalfTmpBuffPtr(0));
              } else { //Perform update instead of write, see Accumulator tyoe
                // Invert Z highpass filtering
                SeparableUpsampledConvolutionEngine3D<T,
                    SubsampledAccumulatorUpdate,
                    typename DTWaveletSchemeT::i_hni
                  >::PerformUpsampledFilteringZRef(
                  this->m_coeff->GetScaleShape(l-1).at(0),
                  this->m_coeff->GetScaleShape(l-1).at(1),
                  this->m_coeff->GetScaleShape(l).at(2),
                  this->m_coeff->GetScaleShape(l-1).at(2),
                  outlow,
                  this->m_coeff->GetHalfTmpBuffPtr(0));
              }
            }
          }
          std::swap(inlow,outlow);
        }
      }

      //Now reconstructing level 1
      int l=1;
      // last filtering stage: merge all 8 octants
      for (int zImagStatus : {0,1} ) {
        for (int zHighStatus : {0,1} ) {
          for (int yImagStatus : {0,1} ) {
            for (int yHighStatus : {0,1} ) {

	      /********************
	       * X Filtering step *
	       ********************/
              // Filtering 4 buffers at once
	      int bandIdxReal = 4*zImagStatus+2*yImagStatus;
	      int bandIdxImag = bandIdxReal+1;
              auto sBandCalc = [=](auto xIdx, auto bandIdx){
                return this->m_coeff->GetHighSubspacePtr(
                  l-1 ,4*zHighStatus+2*yHighStatus+xIdx-1, bandIdx);
              };
	      T *lowReal, *lowImag;
              if ((zHighStatus==0)&&(yHighStatus==0)) {
                if (l==this->m_level) {
                  lowReal = this->m_coeff->GetLowSubspacePtr(l-1,bandIdxReal);
                  lowImag = this->m_coeff->GetLowSubspacePtr(l-1,bandIdxImag);
                } else {
                  lowReal = this->m_coeff->GetHalfTmpBuffPtr(2);
                  lowImag = this->m_coeff->GetHalfTmpBuffPtr(2);
                }
              } else {
                lowReal = sBandCalc(0, bandIdxReal);
                lowImag = sBandCalc(0, bandIdxImag);
              }
    	      // Merge Band 4z+2y+ 0 and 1 in the X direction
	      SeparableUpsampledConvolutionEngine3D<T,
                  SubsampledAccumulator,
		  typename DTWaveletSchemeT::i_l0r,
		  typename DTWaveletSchemeT::i_l0i,
		  typename DTWaveletSchemeT::i_h0r,
		  typename DTWaveletSchemeT::i_h0i
		>::PerformUpsampledFilteringXRef(
		  this->m_coeff->GetScaleShape(l).at(0),
		  this->m_coeff->GetScaleShape(l-1).at(0),
		  this->m_coeff->GetScaleShape(l).at(1),
		  this->m_coeff->GetScaleShape(l).at(2),
		  this->m_coeff->GetHalfTmpBuffPtr(1),
		  lowReal, //in
		  lowImag, //in
                  sBandCalc(1, bandIdxReal), //in
                  sBandCalc(1, bandIdxImag)); //in


	      /********************
	       * Y Filtering step *
	       ********************/
              if (yImagStatus==0) {
		if (yHighStatus==0) {
		  // Invert Y lowpass filtering only
		  SeparableUpsampledConvolutionEngine3D<T,
		      SubsampledAccumulator,
		      typename DTWaveletSchemeT::i_l0r
		    >::PerformUpsampledFilteringYRef(
		    this->m_coeff->GetScaleShape(l-1).at(0),
		    this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l-1).at(1),
		    this->m_coeff->GetScaleShape(l).at(2),
		    this->m_coeff->GetHalfTmpBuffPtr(0),
		    this->m_coeff->GetHalfTmpBuffPtr(1)); 
		} else { //Perform update instead of write, see Accumulator type
		  // Invert Y highpass filtering only
		  SeparableUpsampledConvolutionEngine3D<T,
		      SubsampledAccumulatorUpdate,
		      typename DTWaveletSchemeT::i_h0r
		    >::PerformUpsampledFilteringYRef(
		    this->m_coeff->GetScaleShape(l-1).at(0),
		    this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l-1).at(1),
		    this->m_coeff->GetScaleShape(l).at(2),
		    this->m_coeff->GetHalfTmpBuffPtr(0),
		    this->m_coeff->GetHalfTmpBuffPtr(1));
		}
              } else {
		if (yHighStatus==0) {
		  // Invert Y lowpass filtering only
		  SeparableUpsampledConvolutionEngine3D<T,
		      SubsampledAccumulatorUpdate,
		      typename DTWaveletSchemeT::i_l0i
		    >::PerformUpsampledFilteringYRef(
		    this->m_coeff->GetScaleShape(l-1).at(0),
		    this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l-1).at(1),
		    this->m_coeff->GetScaleShape(l).at(2),
		    this->m_coeff->GetHalfTmpBuffPtr(0),
		    this->m_coeff->GetHalfTmpBuffPtr(1)); 
		} else { //Perform update instead of write, see Accumulator type
		  // Invert Y highpass filtering only
		  SeparableUpsampledConvolutionEngine3D<T,
		      SubsampledAccumulatorUpdate,
		      typename DTWaveletSchemeT::i_h0i
		    >::PerformUpsampledFilteringYRef(
		    this->m_coeff->GetScaleShape(l-1).at(0),
		    this->m_coeff->GetScaleShape(l).at(1),
		    this->m_coeff->GetScaleShape(l-1).at(1),
		    this->m_coeff->GetScaleShape(l).at(2),
		    this->m_coeff->GetHalfTmpBuffPtr(0),
		    this->m_coeff->GetHalfTmpBuffPtr(1));
		}
              }
            }
          }

	  /********************
	   * Z Filtering step *
	   ********************/
	  if (zImagStatus==0) {
	    if (zHighStatus==0) {
	      // Invert Y lowpass filtering only
	      SeparableUpsampledConvolutionEngine3D<T,
		  SubsampledAccumulator,
		  typename DTWaveletSchemeT::i_l0r
		>::PerformUpsampledFilteringZRef(
		this->m_coeff->GetScaleShape(l-1).at(0),
		this->m_coeff->GetScaleShape(l-1).at(1),
		this->m_coeff->GetScaleShape(l).at(2),
		this->m_coeff->GetScaleShape(l-1).at(2),
		this->m_image,
		this->m_coeff->GetHalfTmpBuffPtr(0)); 
	    } else { //Perform update instead of write, see Accumulator type
	      // Invert Y highpass filtering only
	      SeparableUpsampledConvolutionEngine3D<T,
		  SubsampledAccumulatorUpdate,
		  typename DTWaveletSchemeT::i_h0r
		>::PerformUpsampledFilteringZRef(
		this->m_coeff->GetScaleShape(l-1).at(0),
		this->m_coeff->GetScaleShape(l-1).at(1),
		this->m_coeff->GetScaleShape(l).at(2),
		this->m_coeff->GetScaleShape(l-1).at(2),
		this->m_image,
		this->m_coeff->GetHalfTmpBuffPtr(0));
	    }
	  } else {
	    if (zHighStatus==0) {
	      // Invert Y lowpass filtering only
	      SeparableUpsampledConvolutionEngine3D<T,
		  SubsampledAccumulatorUpdate,
		  typename DTWaveletSchemeT::i_l0i
		>::PerformUpsampledFilteringZRef(
		this->m_coeff->GetScaleShape(l-1).at(0),
		this->m_coeff->GetScaleShape(l-1).at(1),
		this->m_coeff->GetScaleShape(l).at(2),
		this->m_coeff->GetScaleShape(l-1).at(2),
		this->m_image,
		this->m_coeff->GetHalfTmpBuffPtr(0)); 
	    } else { //Perform update instead of write, see Accumulator type
	      // Invert Y highpass filtering only
	      SeparableUpsampledConvolutionEngine3D<T,
		  SubsampledAccumulatorUpdate,
		  typename DTWaveletSchemeT::i_h0i
		>::PerformUpsampledFilteringZRef(
		this->m_coeff->GetScaleShape(l-1).at(0),
		this->m_coeff->GetScaleShape(l-1).at(1),
		this->m_coeff->GetScaleShape(l).at(2),
		this->m_coeff->GetScaleShape(l-1).at(2),
		this->m_image,
		this->m_coeff->GetHalfTmpBuffPtr(0));
	    }
	  }
        }
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
using PackedDTContainer3D =
  DTCoeffContainer3D<T,std::vector<T>>;

// Aliasing ugly types into more simple ones
template<typename T>
using dtwAnto97QSHIFT6_3D = 
  DTWavelet3D<T,PackedDTContainer3D<T>,dtwAnto97QSHIFT6<T>>;

/** \struct DB3DWt
 * \brief Utility struct that allow to instanciate all 3D wavelets at once
 *
 * \author Thibault Notargiacomo
 */
template<typename T>
struct DB3DWt {
 Daub2_3D<T> daub2_3D;
 Daub3_3D<T> daub3_3D;
 Daub4_3D<T> daub4_3D;
 Daub5_3D<T> daub5_3D;
 Anto97_BiOrth_3D<T> anto97_BiOrth_3D;
 QSHIFT6_Orth_3D<T> QShift6_Orth_3D;
 REVERSE_QSHIFT6_Orth_3D<T> Reverse_Qshift6_Orth_3D;
 dtwAnto97QSHIFT6_3D<T> dtwAnto97QShift6_3D; 
};



#endif //WAVELET3D_H
