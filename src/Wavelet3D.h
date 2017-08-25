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
    T* outlowX;
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
          
          // If this is the low freq projection
          if ((zFiltIdx==0)&&(yFiltIdx==0)) {
			if (l+1==this->m_level) {
			  outlowX=this->m_coeff->GetLowSubspacePtr(l);
			} else {
			  outlowX=this->m_coeff->GetOutLowTmpBuffPtr();
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
      inlow = outlowX;
    }
    return 1;
  }


  /// Backward wavelet transform: transpose of the forward transpose
  virtual int backward() {
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
            inlowX=this->m_coeff->GetLowSubspacePtr(l-1);
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

		//Update output buffer destination TODO TN URGENT:
        //strategy does not work because one overwrites data for l>1
		T* outlow;
		if (l<=1) {
		  outlow=this->m_image;
		} else {
		  outlow=this->m_coeff->GetLowSubspacePtr(l-2);
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
			  outlow,
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
			  outlow,
			  this->m_coeff->GetHalfTmpBuffPtr(0));
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
