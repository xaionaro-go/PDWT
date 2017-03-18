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
  Filter()=default;

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
template<typename T, int TAP_SIZE_LEFT, int TAP_SIZE_RIGHT>
struct Filter : public GenericFilter<T,TAP_SIZE_LEFT,TAP_SIZE_RIGHT> {
  /// Defaulted constructor
  Filter()=default;
  
  /// Actual storage for the filter
  static const T Buf[Filter<T,TAP_SIZE_LEFT,TAP_SIZE_RIGHT>::TapSize];
};

struct wfilter {
  char wname[16];
  int hlen;
  DTYPE* f_l; // Forward lowpass
  DTYPE* f_h;  // Forward highpass
  DTYPE* i_l; // Inverse lowpass
  DTYPE* i_h; // Inverse highpass
};

wfilter all_filters[72];

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
