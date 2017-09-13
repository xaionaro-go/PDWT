#ifndef WAVELET3D_H
#define WAVELET3D_H

// Local
#include "Wavelet.h"

// STL
#include <vector>

// Local
#include "coeffContainer.h"
#include "filters.h"
#include "separable.h"
#include "vectorization/vectorization.h"


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
    T* outlowZReallowYReallowXReal;
    T* outlowZReallowYReallowXImag;
    T* outlowZReallowYImaglowXReal;
    T* outlowZReallowYImaglowXImag;
    T* outlowZImaglowYReallowXReal;
    T* outlowZImaglowYReallowXImag;
    T* outlowZImaglowYImaglowXReal;
    T* outlowZImaglowYImaglowXImag;
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

						} else { 
              //Z Real Filters high
						}
          } else {
            //Z Imag Filters
						if (zFiltIdx==0) {
              //Z Imag Filters low

						} else {
              //Z Imag Filters high

						}
          }
   
					/********************
					 * Y Filtering step *
					 ********************/
					//Filter with first yReal filter, and then yImag filter
					for (int yImagStatus : {0,1}) {

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
// dtwAnto97QSHIFT6_3D<T> dtwAnto97QShift6_3D; 
};



#endif //WAVELET3D_H
