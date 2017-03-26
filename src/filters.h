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
#include <string>

// Local
#include "vectorization.h"

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
  
  /// Total size of the filter, in number of elements
  constexpr static int TapSize =
    TAP_SIZE_LEFT + TAP_SIZE_RIGHT + 1; //+1 = the center pixel
  
  /// Redefine template parameters
  constexpr static int VecSize = sizeof(VectorType)/sizeof(T);
  
  /// How many vector are needed to load a single filter support
  constexpr static int NbVecPerFilt = (TapSize+VecSize-1)/(VecSize);
  
  /// Compile time definition of filter feature: left tap number
  constexpr static int TapSizeLeft = TAP_SIZE_LEFT;
  
  /// Compile time definition of filter feature: right tap number
  constexpr static int TapSizeRight = TAP_SIZE_RIGHT;
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
  /// Defaulted constructor
  Filter()=default;
  
  // Actual storage for the filter
  static constexpr T Buf[GenericFilter<
    T,TAP_SIZE_LEFT,TAP_SIZE_RIGHT>::TapSize]={0};
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
};

//Semi specialization, type agnostic
template<typename T>
struct Filter<T,1,2,filterDB,filterDB::DB2_L> {
  static constexpr T Buff[4] = { -0.12940952255092145, 0.22414386804185735,
    0.836516303737469, 0.48296291314469025 };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,1,2,filterDB,filterDB::DB2_H> {
  static constexpr T Buff[4] = { -0.48296291314469025, 0.836516303737469,
    -0.22414386804185735, -0.12940952255092145 };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,2,1,filterDB,filterDB::DB2_I_L> {
  static constexpr T Buff[4] = { 0.48296291314469025, 0.836516303737469,
    0.22414386804185735, -0.12940952255092145 };
};
//Semi specialization, type agnostic
template<typename T>
struct Filter<T,2,1,filterDB,filterDB::DB2_I_H> {
  static constexpr T Buff[4] = { -0.12940952255092145, -0.22414386804185735,
    0.836516303737469, -0.48296291314469025 };
};

/// The Daubechies2 wavelet system
template<typename T>
using Daub2<T> = wFilter<Filter<T,1,2,filterDB,filterDB::DB2_L>,
    Filter<T,1,2,filterDB,filterDB::DB2_H>,
    Filter<T,2,1,filterDB,filterDB::DB2_I_L>,
    Filter<T,2,1,filterDB,filterDB::DB2_I_H>>;


/*
DB2_L[4];
DB2_H[4];
DB2_I_L[4];
DB2_I_H[4];
DB3_L[6];
DB3_H[6];
DB3_I_L[6];
DB3_I_H[6];
DB4_L[8];
DB4_H[8];
DB4_I_L[8];
DB4_I_H[8];
DB5_L[10];
DB5_H[10];
DB5_I_L[10];
DB5_I_H[10];
DB6_L[12];
DB6_H[12];
DB6_I_L[12];
DB6_I_H[12];
DB7_L[14];
DB7_I_L[14];
DB7_H[14];
DB7_I_H[14];
DB8_L[16];
DB8_I_L[16];
DB8_H[16];
DB8_I_H[16];
DB9_L[18];
DB9_I_L[18];
DB9_H[18];
DB9_I_H[18];
DB10_L[20];
DB10_I_L[20];
DB10_H[20];
DB10_I_H[20];
DB11_L[22];
DB11_I_L[22];
DB11_H[22];
DB11_I_H[22];
DB12_L[24];
DB12_I_L[24];
DB12_H[24];
DB12_I_H[24];
DB13_L[26];
DB13_I_L[26];
DB13_H[26];
DB13_I_H[26];
DB14_L[28];
DB14_I_L[28];
DB14_H[28];
DB14_I_H[28];
DB15_L[30];
DB15_I_L[30];
DB15_H[30];
DB15_I_H[30];
DB16_L[32];
DB16_I_L[32];
DB16_H[32];
DB16_I_H[32];
DB17_L[34];
DB17_I_L[34];
DB17_H[34];
DB17_I_H[34];
DB18_L[36];
DB18_I_L[36];
DB18_H[36];
DB18_I_H[36];
DB19_L[38];
DB19_I_L[38];
DB19_H[38];
DB19_I_H[38];
DB20_L[40];
DB20_I_L[40];
DB20_H[40];
DB20_I_H[40];
SYM2_L[4];
SYM2_I_L[4];
SYM2_H[4];
SYM2_I_H[4];
SYM3_L[6];
SYM3_I_L[6];
SYM3_H[6];
SYM3_I_H[6];
SYM4_L[8];
SYM4_I_L[8];
SYM4_H[8];
SYM4_I_H[8];
SYM5_L[10];
SYM5_I_L[10];
SYM5_H[10];
SYM5_I_H[10];
SYM6_L[12];
SYM6_I_L[12];
SYM6_H[12];
SYM6_I_H[12];
SYM7_L[14];
SYM7_I_L[14];
SYM7_H[14];
SYM7_I_H[14];
SYM8_L[16];
SYM8_I_L[16];
SYM8_H[16];
SYM8_I_H[16];
SYM9_L[18];
SYM9_I_L[18];
SYM9_H[18];
SYM9_I_H[18];
SYM10_L[20];
SYM10_I_L[20];
SYM10_H[20];
SYM10_I_H[20];
SYM11_L[22];
SYM11_I_L[22];
SYM11_H[22];
SYM11_I_H[22];
SYM12_L[24];
SYM12_I_L[24];
SYM12_H[24];
SYM12_I_H[24];
SYM13_L[26];
SYM13_I_L[26];
SYM13_H[26];
SYM13_I_H[26];
SYM14_L[28];
SYM14_I_L[28];
SYM14_H[28];
SYM14_I_H[28];
SYM15_L[30];
SYM15_I_L[30];
SYM15_H[30];
SYM15_I_H[30];
SYM16_L[32];
SYM16_I_L[32];
SYM16_H[32];
SYM16_I_H[32];
SYM17_L[34];
SYM17_I_L[34];
SYM17_H[34];
SYM17_I_H[34];
SYM18_L[36];
SYM18_I_L[36];
SYM18_H[36];
SYM18_I_H[36];
SYM19_L[38];
SYM19_I_L[38];
SYM19_H[38];
SYM19_I_H[38];
SYM20_L[40];
SYM20_I_L[40];
SYM20_H[40];
SYM20_I_H[40];
COIF1_L[6];
COIF1_I_L[6];
COIF1_H[6];
COIF1_I_H[6];
COIF2_L[12];
COIF2_I_L[12];
COIF2_H[12];
COIF2_I_H[12];
COIF3_L[18];
COIF3_I_L[18];
COIF3_H[18];
COIF3_I_H[18];
COIF4_L[24];
COIF4_I_L[24];
COIF4_H[24];
COIF4_I_H[24];
COIF5_L[30];
COIF5_I_L[30];
COIF5_H[30];
COIF5_I_H[30];
BIOR1_3_L[6];
BIOR1_3_I_L[6];
BIOR1_3_H[6];
BIOR1_3_I_H[6];
BIOR1_5_L[10];
BIOR1_5_I_L[10];
BIOR1_5_H[10];
BIOR1_5_I_H[10];
BIOR2_2_L[6];
BIOR2_2_I_L[6];
BIOR2_2_H[6];
BIOR2_2_I_H[6];
BIOR2_4_L[10];
BIOR2_4_I_L[10];
BIOR2_4_H[10];
BIOR2_4_I_H[10];
BIOR2_6_L[14];
BIOR2_6_I_L[14];
BIOR2_6_H[14];
BIOR2_6_I_H[14];
BIOR2_8_L[18];
BIOR2_8_I_L[18];
BIOR2_8_H[18];
BIOR2_8_I_H[18];
BIOR3_1_L[4];
BIOR3_1_I_L[4];
BIOR3_1_H[4];
BIOR3_1_I_H[4];
BIOR3_3_L[8];
BIOR3_3_I_L[8];
BIOR3_3_H[8];
BIOR3_3_I_H[8];
BIOR3_5_L[12];
BIOR3_5_I_L[12];
BIOR3_5_H[12];
BIOR3_5_I_H[12];
BIOR3_7_L[16];
BIOR3_7_I_L[16];
BIOR3_7_H[16];
BIOR3_7_I_H[16];
BIOR3_9_L[20];
BIOR3_9_I_L[20];
BIOR3_9_H[20];
BIOR3_9_I_H[20];
BIOR4_4_L[10];
BIOR4_4_I_L[10];
BIOR4_4_H[10];
BIOR4_4_I_H[10];
BIOR5_5_L[12];
BIOR5_5_H[12];
BIOR5_5_I_L[12];
BIOR5_5_I_H[12];
BIOR6_8_L[18];
BIOR6_8_H[18];
BIOR6_8_I_L[18];
BIOR6_8_I_H[18];
RBIOR1_3_L[6];
RBIOR1_3_I_L[6];
RBIOR1_3_H[6];
RBIOR1_3_I_H[6];
RBIOR1_5_L[10];
RBIOR1_5_I_L[10];
RBIOR1_5_H[10];
RBIOR1_5_I_H[10];
RBIOR2_2_L[6];
RBIOR2_2_I_L[6];
RBIOR2_2_H[6];
RBIOR2_2_I_H[6];
RBIOR2_4_L[10];
RBIOR2_4_I_L[10];
RBIOR2_4_H[10];
RBIOR2_4_I_H[10];
RBIOR2_6_L[14];
RBIOR2_6_I_L[14];
RBIOR2_6_H[14];
RBIOR2_6_I_H[14];
RBIOR2_8_L[18];
RBIOR2_8_I_L[18];
RBIOR2_8_H[18];
RBIOR2_8_I_H[18];
RBIOR3_1_L[4];
RBIOR3_1_I_L[4];
RBIOR3_1_H[4];
RBIOR3_1_I_H[4];
RBIOR3_3_L[8];
RBIOR3_3_I_L[8];
RBIOR3_3_H[8];
RBIOR3_3_I_H[8];
RBIOR3_5_L[12];
RBIOR3_5_I_L[12];
RBIOR3_5_H[12];
RBIOR3_5_I_H[12];
RBIOR3_7_L[16];
RBIOR3_7_I_L[16];
RBIOR3_7_H[16];
RBIOR3_7_I_H[16];
RBIOR3_9_L[20];
RBIOR3_9_I_L[20];
RBIOR3_9_H[20];
RBIOR3_9_I_H[20];
RBIOR4_4_L[10];
RBIOR4_4_I_L[10];
RBIOR4_4_H[10];
RBIOR4_4_I_H[10];
RBIOR5_5_L[12];
RBIOR5_5_I_L[12];
RBIOR5_5_H[12];
RBIOR5_5_I_H[12];
RBIOR6_8_L[18];
RBIOR6_8_I_L[18];
RBIOR6_8_H[18];
RBIOR6_8_I_H[18];
HAAR_L[4];
HAAR_H[4];
HAAR_I_L[4];
HAAR_I_H[4];
*/

#endif
