///
/// Wavelets coefficients, extracted from http://wavelets.pybytes.com
///

#ifndef FILTERS_H
#define FILTERS_H

/*
TODO TN: a specialiser et supprimer
#define DTYPE float
#define cublas_asum cublasSasum
#define cublas_nrm2 cublasSnrm2
#define cublas_scal cublasSscal
#define cublas_axpy cublasSaxpy
#else
//~ typedef double DTYPE;
#define DTYPE double
#define cublas_asum cublasDasum
#define cublas_nrm2 cublasDnrm2
#define cublas_scal cublasDscal
#define cublas_axpy cublasDaxpy
#endif*/


// STL
#include <array>
#include <string>

// Local
#include "vectorization/vectorization.h"

/** \struct GenericFilter
 * \brief Generic interface that should be implemented by filters. It has been
 * designed to provide various metaprogramming helper in order to setup
 * explicit vectorization or enable advanced compiler inlining/optimization
 * 
 * TAP_SIZE_LEFT does only account for number of elements at left,
 * without the current one.
 * TAP_SIZE_RIGHT does only account for number of elements at right,
 * without the current one.
 *
 * \author Thibault Notargiacomo
 */
template<typename T, int TAP_SIZE_LEFT, int TAP_SIZE_RIGHT>
struct GenericFilter {
  /// Defaulted constructor
  GenericFilter()=default;

  /// Typedef for scalar main type
  using ScalarType=T;
  
  /// Typedef for vectorization type
  using VectorType=PackType<T>;
  
  /// Redefine template parameters
  constexpr static int VecSize = sizeof(VectorType)/sizeof(T);

  /// Compile time definition of filter feature: left tap number
  constexpr static int TapSizeLeft = TAP_SIZE_LEFT;
  
  /// Compile time definition of filter feature: right tap number
  constexpr static int TapSizeRight = TAP_SIZE_RIGHT;

  /// Total size of the filter, in number of elements
  constexpr static int TapSize =
    TapSizeLeft + TapSizeRight + 1; //+1 = the center pixel

  /// For specific cases in conv with subsampled data
  constexpr static int TapHalfSize = TapSize/2;

  /// For specific cases in conv with subsampled data
  constexpr static bool IsHalfSizeOdd = TapHalfSize&1;

  /// for conv with subsampled data: half tap number
  constexpr static int TapHalfSizeLeft = TapHalfSize/2;
  
  /// for conv with subsampled data: half tap number
  constexpr static int TapHalfSizeRight = TapHalfSizeLeft-(IsHalfSizeOdd?0:1);

  /// How many vector are needed to load a single filter support
  constexpr static int NbVecPerFilt = (TapSize+VecSize-1)/(VecSize);
 };

/** \struct Filter
 * \bried We define an inheritance of the fully generic filter for a structure
 * that will be able to handle its own storage for the filter
 *
 * \author Thibault Notargiacomo
 */
template<typename T, int TAP_SIZE_LEFT, int TAP_SIZE_RIGHT,
  typename FilterEnum, FilterEnum FilterVal>
struct Filter : public GenericFilter<T,TAP_SIZE_LEFT,TAP_SIZE_RIGHT> {
  using ParenT =  GenericFilter<T,TAP_SIZE_LEFT,TAP_SIZE_RIGHT>; 

  /// Defaulted constructor
  Filter()=default;
  
  // Actual storage for the filter
  static const constexpr std::array<T,ParenT::TapSize> Buff = {0};
};

/** \struct wFilter
 * \brief A group of filters that defines a full forward and backward filtering
 * process
 *
 * \author Pierre Paleo and Thibault Notargiacomo
 */
template<class ForwardLowT, class ForwardHighT, class InverseLowT,
  class InverseHighT>
struct wFilter {
  /// A small string that defines the wavelet system name
  std::string wname;
  
  /// Forward lowpass filter
  using f_l = ForwardLowT;
  /// Forward highpasas filter
  using f_h = ForwardHighT;    
  /// Inverse lowpass filter
  using i_l = InverseLowT; 
  /// Inverse highpass filter
  using i_h = InverseHighT;
};

/** \struct cwFilter
 * \brief A group of filters that defines a full forward and backward complex
 * wavelet filter process
 *
 * \author Thibault Notargiacomo
 */
template<class ForwardLowStage0RealT,
         class ForwardLowStage0ImagT,
         class ForwardHighStage0RealT,
         class ForwardHighStage0ImagT,
         class ForwardLowStageNRealT,
         class ForwardLowStageNImagT,
         class ForwardHighStageNRealT,
         class ForwardHighStageNImagT,
         class InverseLowStage0RealT,
         class InverseLowStage0ImagT,
         class InverseHighStage0RealT,
         class InverseHighStage0ImagT,
         class InverseLowStageNRealT,
         class InverseLowStageNImagT,
         class InverseHighStageNRealT,
         class InverseHighStageNImagT>
struct cwFilter {
  /// A small string that defines the wavelet system name
  std::string wname;
  
  /// Forward lowpass filter stage 0 real part
  using f_l0r = ForwardLowStage0RealT;
  /// Forward lowpass filter stage 0 imaginary part
  using f_l0i = ForwardLowStage0ImagT;
  /// Forward highpass filter stage 0 real part
  using f_h0r = ForwardHighStage0RealT;
  /// Forward highpass filter stage 0 imaginary part
  using f_h0i = ForwardHighStage0ImagT;
  /// Forward lowpass filter stage n real part
  using f_lnr = ForwardLowStageNRealT;
  /// Forward lowpass filter stage n imaginary part
  using f_lni = ForwardLowStageNImagT;
  /// Forward highpass filter stage n real part
  using f_hnr = ForwardHighStageNRealT;
  /// Forward highpass filter stage n imaginary part
  using f_hni = ForwardHighStageNImagT;
  /// Inverse lowpass filter stage 0 real part
  using i_l0r = InverseLowStage0RealT;
  /// Inverse lowpass filter stage 0 imaginary part
  using i_l0i = InverseLowStage0ImagT;
  /// Inverse highpass filter stage 0 real part
  using i_h0r = InverseHighStage0RealT;
  /// Inverse highpass filter stage 0 imaginary part
  using i_h0i = InverseHighStage0ImagT;
  /// Inverse lowpass filter stage n real part
  using i_lnr = InverseLowStageNRealT;
  /// Inverse lowpass filter stage n imaginary part
  using i_lni = InverseLowStageNImagT;
  /// Inverse highpass filter stage n real part
  using i_hnr = InverseHighStageNRealT;
  /// Inverse highpass filter stage n imaginary part
  using i_hni = InverseHighStageNImagT;
};

enum class filterDB {
  DB2_L,
  DB2_H,
  DB2_I_L,
  DB2_I_H,
  DB3_L,
  DB3_H,
  DB3_I_L,
  DB3_I_H,
  DB4_L,
  DB4_H,
  DB4_I_L,
  DB4_I_H,
  DB5_L,
  DB5_H,
  DB5_I_L,
  DB5_I_H,
  ANTO9_L,
  ANTO7_H,
  ANTO7_I_L,
  ANTO9_I_H
};

//Semi specialization, type agnostic
template<typename T>
struct Filter<T,1,2,filterDB,filterDB::DB2_L> : public
    GenericFilter<T,1,2> {
  static const constexpr std::array<T,4> Buff = {
    -0.12940952255092145,
    0.22414386804185735,
    0.836516303737469,
    0.48296291314469025
  };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,1,2,filterDB,filterDB::DB2_H> : public
    GenericFilter<T,1,2> {
  static const constexpr std::array<T,4> Buff = {
    -0.48296291314469025,
    0.836516303737469,
    -0.22414386804185735,
    -0.12940952255092145
  };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,2,1,filterDB,filterDB::DB2_I_L> : public
    GenericFilter<T,2,1> {
  static const constexpr std::array<T,4> Buff = {
    0.48296291314469025,
    0.836516303737469,
    0.22414386804185735,
    -0.12940952255092145
  };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,2,1,filterDB,filterDB::DB2_I_H> : public
    GenericFilter<T,2,1> {
  static const constexpr std::array<T,4> Buff = {
    -0.12940952255092145,
    -0.22414386804185735,
    0.836516303737469,
    -0.48296291314469025
  };
};
/// The Daubechies2 wavelet system, type agnostic
template<typename T>
using Daub2 = wFilter<
    Filter<T,1,2,filterDB,filterDB::DB2_L>,
    Filter<T,1,2,filterDB,filterDB::DB2_H>,
    Filter<T,2,1,filterDB,filterDB::DB2_I_L>,
    Filter<T,2,1,filterDB,filterDB::DB2_I_H>>;


//Semi specialization, type agnostic
template<typename T>
struct Filter<T,2,3,filterDB,filterDB::DB3_L> : public
    GenericFilter<T,2,3> {
  static const constexpr std::array<T,6> Buff = {
    0.035226291882100656,
    -0.08544127388224149,
    -0.13501102001039084,
    0.4598775021193313, 
    0.8068915093133388,
    0.3326705529509569
  };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,2,3,filterDB,filterDB::DB3_H> : public
    GenericFilter<T,2,3> {
  static const constexpr std::array<T,6> Buff = {
    -0.3326705529509569,
    0.8068915093133388,
    -0.4598775021193313,
    -0.13501102001039084,
    0.08544127388224149,
    0.035226291882100656
  };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,3,2,filterDB,filterDB::DB3_I_L> : public
    GenericFilter<T,3,2> {
  static const constexpr std::array<T,6> Buff = {
    0.3326705529509569, 
    0.8068915093133388,
    0.4598775021193313,
    -0.13501102001039084,
    -0.08544127388224149,
    0.035226291882100656
  };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,3,2,filterDB,filterDB::DB3_I_H> : public
    GenericFilter<T,3,2> {
  static const constexpr std::array<T,6> Buff = {
    0.035226291882100656,
    0.08544127388224149,
    -0.13501102001039084,
    -0.4598775021193313,
    0.8068915093133388,
    -0.3326705529509569
  };
};
/// The Daubechies3 wavelet system, type agnostic
template<typename T>
using Daub3 = wFilter<
    Filter<T,2,3,filterDB,filterDB::DB3_L>,
    Filter<T,2,3,filterDB,filterDB::DB3_H>,
    Filter<T,3,2,filterDB,filterDB::DB3_I_L>,
    Filter<T,3,2,filterDB,filterDB::DB3_I_H>>;


//Semi specialization, type agnostic
template<typename T>
struct Filter<T,3,4,filterDB,filterDB::DB4_L> : public
    GenericFilter<T,3,4> {
  static const constexpr std::array<T,8> Buff = {
    -0.010597401784997278,
    0.032883011666982945 ,
    0.030841381835986965 ,
    -0.18703481171888114 ,
    -0.02798376941698385 ,
    0.6308807679295904   ,
    0.7148465705525415   ,
    0.23037781330885523
  };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,3,4,filterDB,filterDB::DB4_H> : public
    GenericFilter<T,3,4> {
  static const constexpr std::array<T,8> Buff = {
    -0.23037781330885523 ,
    0.7148465705525415   ,
    -0.6308807679295904  ,
    -0.02798376941698385 ,
    0.18703481171888114  ,
    0.030841381835986965 ,
    -0.032883011666982945,
    -0.010597401784997278
  };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,4,3,filterDB,filterDB::DB4_I_L> : public
    GenericFilter<T,4,3> {
  static const constexpr std::array<T,8> Buff = {
    0.23037781330885523  ,
    0.7148465705525415   ,
    0.6308807679295904   ,
    -0.02798376941698385 ,
    -0.18703481171888114 ,
    0.030841381835986965 ,
    0.032883011666982945 ,
    -0.010597401784997278
  };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,4,3,filterDB,filterDB::DB4_I_H> : public
    GenericFilter<T,4,3> {
  static const constexpr std::array<T,8> Buff = {
    -0.010597401784997278,
    -0.032883011666982945,
    0.030841381835986965 ,
    0.18703481171888114  ,
    -0.02798376941698385 ,
    -0.6308807679295904  ,
    0.7148465705525415   ,
    -0.23037781330885523
  };
};
/// The Daubechies4 wavelet system, type agnostic
template<typename T>
using Daub4 = wFilter<
    Filter<T,3,4,filterDB,filterDB::DB4_L>,
    Filter<T,3,4,filterDB,filterDB::DB4_H>,
    Filter<T,4,3,filterDB,filterDB::DB4_I_L>,
    Filter<T,4,3,filterDB,filterDB::DB4_I_H>>;


//Semi specialization, type agnostic
template<typename T>
struct Filter<T,4,5,filterDB,filterDB::DB5_L> : public
    GenericFilter<T,4,5> {
  static const constexpr std::array<T,10> Buff = {
    0.003335725285001549,
    -0.012580751999015526,
    -0.006241490213011705,
    0.07757149384006515,
    -0.03224486958502952,
    -0.24229488706619015,
    0.13842814590110342,
    0.7243085284385744,
    0.6038292697974729,
    0.160102397974125
 };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,4,5,filterDB,filterDB::DB5_H> : public
    GenericFilter<T,4,5> {
  static const constexpr std::array<T,10> Buff = {
    -0.160102397974125,
    0.6038292697974729,
    -0.7243085284385744,
    0.13842814590110342,
    0.24229488706619015,
    -0.03224486958502952,
    -0.07757149384006515,
    -0.006241490213011705,
    0.012580751999015526,
    0.003335725285001549
 };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,5,4,filterDB,filterDB::DB5_I_L> : public
    GenericFilter<T,5,4> {
  static const constexpr std::array<T,10> Buff = {
    0.160102397974125,
    0.6038292697974729,
    0.7243085284385744,
    0.13842814590110342,
    -0.24229488706619015,
    -0.03224486958502952,
    0.07757149384006515,
    -0.006241490213011705,
    -0.012580751999015526,
    0.003335725285001549
 };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,5,4,filterDB,filterDB::DB5_I_H> : public
    GenericFilter<T,5,4> {
  static const constexpr std::array<T,10> Buff = {
    0.003335725285001549,
    0.012580751999015526,
    -0.006241490213011705,
    -0.07757149384006515,
    -0.03224486958502952,
    0.24229488706619015,
    0.13842814590110342,
    -0.7243085284385744,
    0.6038292697974729,
    -0.160102397974125
 };
};
/// The Daubechies5 wavelet system, type agnostic
template<typename T>
using Daub5 = wFilter<
    Filter<T,4,5,filterDB,filterDB::DB5_L>,
    Filter<T,4,5,filterDB,filterDB::DB5_H>,
    Filter<T,5,4,filterDB,filterDB::DB5_I_L>,
    Filter<T,5,4,filterDB,filterDB::DB5_I_H>>;


//Semi specialization, type agnostic
template<typename T>
struct Filter<T,5,6,filterDB,filterDB::ANTO9_L> : public
    GenericFilter<T,5,6> {
  static const constexpr std::array<T,12> Buff = {
    0.,
    0.02674875741081,
    -0.01686411844287,
    -0.07822326652899,
    0.26686411844288,
    0.60294901823636,
    0.26686411844287,
    -0.07822326652899,
    -0.01686411844287,
    0.02674875741081,
    0.,0.
 };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,5,6,filterDB,filterDB::ANTO7_H> : public
    GenericFilter<T,5,6> {
  static const constexpr std::array<T,12> Buff = {
    0.,0.,0.,
    0.04563588155712,
    -0.02877176311425,
    -0.29563588155712,
    0.55754352622850,
    -0.29563588155713,
    -0.02877176311425,
    0.04563588155712,
    0.,0.
 };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,6,5,filterDB,filterDB::ANTO7_I_L> : public
    GenericFilter<T,6,5> {
  static const constexpr std::array<T,12> Buff = {
    0.,0.,0.,
    -0.09127176311424,
    -0.05754352622850,
    0.59127176311424,
    1.11508705245700,
    0.59127176311426,
    -0.05754352622850,
    -0.09127176311424,
    0.,0.
 };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,6,5,filterDB,filterDB::ANTO9_I_H> : public
    GenericFilter<T,6,5> {
  static const constexpr std::array<T,12> Buff = {
    0,
    0.05349751482162,
    0.03372823688574,
    -0.15644653305798,
    -0.53372823688576,
    1.20589803647272,
    -0.53372823688574,
    -0.15644653305798,
    0.03372823688574,
    0.05349751482162,
    0,0
 };
};
/// The (9,7) tap bi-orthogonal Antonini filter, type agnostic
template<typename T>
using Anto97_BiOrth = wFilter<
    Filter<T,5,6,filterDB,filterDB::ANTO9_L>,
    Filter<T,5,6,filterDB,filterDB::ANTO7_H>,
    Filter<T,6,5,filterDB,filterDB::ANTO7_I_L>,
    Filter<T,6,5,filterDB,filterDB::ANTO9_I_H>>;
/*    Filter<T,4,4,filterDB,filterDB::ANTO9_L>,
    Filter<T,3,3,filterDB,filterDB::ANTO7_H>,
    Filter<T,3,3,filterDB,filterDB::ANTO7_I_L>,
    Filter<T,4,4,filterDB,filterDB::ANTO9_I_H>>;*/

template <typename T>
const std::array<T,12> Filter<T,5,6,filterDB,filterDB::ANTO9_L>::Buff;
template <typename T>
const std::array<T,12> Filter<T,5,6,filterDB,filterDB::ANTO7_H>::Buff;

/*See ../../Pix3DGIT/tools/Volumix2/src/Cuda/CudaConstantMemoryHelper.cpp 
and ../../Pix3DGIT/tools/Volumix2/src/Cuda/DualTreeComplexWaveletOperator/DTCWT.cu.h
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,2,3,filterDB,filterDB::QSHIFT6_L> : public
    GenericFilter<T,2,3> {
  static const constexpr std::array<T,6> Buff = {
 };
};
/// The 6 tap orthogonal Q-Shift filter, type agnostic
template<typename T>
using QSHIFT6_Orth = wFilter<
    Filter<T,4,4,filterDB,filterDB::QSHIFT6_L>,
    Filter<T,3,3,filterDB,filterDB::QSHIFT6_H>,
    Filter<T,3,3,filterDB,filterDB::QSHIFT6_I_L>,
    Filter<T,3,3,filterDB,filterDB::QSHIFT6_I_H>>;
*/



#endif
