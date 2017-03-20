#ifndef WT_H
#define WT_H

// STL
#include <string>
#include <list>

// Local

/**
 * Possible states of the Wavelet class.
 * It prevents, for example, W.inverse() from being run twice
 * (since W.d_coeffs[0] is modified)
 */
enum class w_state {
  /// The class has just been initialized (coeffs not computed)
  W_INIT,
  /// W.forward() has just been performed (coeffs computed)
  W_FORWARD,
  /// W.inverse() has just been performed (d_image modified, coeffs modified !)
  W_INVERSE,
  /// The coefficients have been modified
  W_THRESHOLD,
  /// Error when creating the Wavelets instance
  W_CREATION_ERROR,
  /// Error when computing the forward transform
  W_FORWARD_ERROR,
  /// Error when computing the inverse transform
  W_INVERSE_ERROR,
  /// Error when thresholding the coefficients
  W_THRESHOLD_ERROR
};

/** \class Wavelets
 * \brief The wavelet class manage the lifecycle of the wavelet transform
 * coefficients.
 *
 * \author Pierre Paleo
 */
template<typename T>
class Wavelets {
 public:
  /// Defaulted constructor
  Wavelets()=default;
  /// Constructor : Wavelets from image
  Wavelets(T* img, int Nr, int Nc, const char* wname, int levels,
      int memisonhost=1, bool do_cycle_spinning=false, int do_swt=0,
      int ndim=2); 
  /// Copy constructor
  Wavelets(const Wavelets &W);
  /// Default destructor
  virtual ~Wavelets();

  /// Forward wavelet tranform
  void forward();
  void soft_threshold(T beta, int do_thresh_appcoeffs = 0, int normalize = 0);
  void hard_threshold(T beta, int do_thresh_appcoeffs = 0, int normalize = 0);
  void shrink(T beta, int do_thresh_appcoeffs = 1);
  void proj_linf(T beta, int do_thresh_appcoeffs = 1);
  void circshift(int sr, int sc, int inplace = 1);
  void inverse();
  T norm2sq();
  T norm1();
  int get_image(T* img);
  void print_informations();
  int get_coeff(T* coeff, int num);
  void set_image(T* img, int mem_is_on_device = 0);
  void set_coeff(T* coeff, int num, int mem_is_on_device = 0);
  int set_filters_forward(char* filtername, uint len, T* filter1, T* filter2, T* filter3 = NULL, T* filter4 = NULL);
  int set_filters_inverse(T* filter1, T* filter2, T* filter3 = NULL, T* filter4 = NULL);

  int add_wavelet(Wavelets W, T alpha=1.0f);
public:
  /// Image (input or result of reconstruction), on device
  T* m_image;
  /// Wavelet coefficients, on device
  std::unique_ptr<CoeffContainerT> m_coeffs;
  /// Temporary device array (to avoid multiple malloc/free)
  //T* d_tmp;
  
  /// Current shift for the cycle spinning process
  std::list<int> m_currentShift;
  
  /// Wavelet name
  std::string m_name;
  /** If cycle spinning is enabled, use image shifting for translation
   * invariant denoising
   */
  bool m_doCycleSpinning;

  /// Informations about the wavelet tranform
  w_info winfos;

  /// Current state of the wavelet tranform
  w_state state;


};


#endif

