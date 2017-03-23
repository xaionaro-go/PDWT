#ifndef WT_H
#define WT_H

// STL
#include <string>
#include <list>

// Local

/** \struct w_info
 * \brief Description of the workload
 *
 * \author Pierre Paleo
 */
struct w_info {
  /// Number of dimensions. For now only 2D and (batched) 1D are supported
  int ndims;
  /// Number of rows of the image (for 1D : Nr = 1)
  int Nr;
  /// Number of columns in the image
  int Nc;
  /// Number of slices in the image (3D)
  int Ns;
  /// Number of decomposition levels
  int nlevels;
  /// Do Stationary (Undecimated) Wavelet Transform
  int do_swt;              
  /// "Filter" length, deprecated in the compile time engine
  int hlen;
};

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
template<typename T, class CoeffContainerT, class WaveletSchemeT>
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
  /// Backward wavelet transform: transpose of the forward transpose
  void backward();
  /// Inverse of the wavelet tranform
  void inverse();
 
  /// Soft thresholding: proximity operator for L1 norm
  //void soft_threshold(T beta, int do_thresh_appcoeffs = 0, int normalize = 0);
  /// Hard thresholding: best k-term approximation for wavelets 
  //void hard_threshold(T beta, int do_thresh_appcoeffs = 0, int normalize = 0);
  /// Wavelet shrinkage
  //void shrink(T beta, int do_thresh_appcoeffs = 1);
  /// Projection on the \f$ l-\infty \f$ ball
  //void proj_linf(T beta, int do_thresh_appcoeffs = 1);
  /// Circular shift
  //void circshift(int sr, int sc, int inplace = 1);
  

  /// Compute the \f$ l-2 \f$ norm of the vector of coefficients
  //T norm2sq();
  /// Compute the \f$ l-1 \f$ norm of the vector of coefficients
  //T norm1();
  /// Get image pointer
  int get_image(T* img);
  /// Print wavelet transform informations
  void print_informations();
  /// Return a pointer to wavelet coefficients
  int get_coeff(T* coeff, int num);
  /// Set input image to be transformed in the wavelet domain
  void set_image(T* img, int mem_is_on_device = 0);
  /// Set coefficients to be reconstructed into an image
  void set_coeff(T* coeff, int num, int mem_is_on_device = 0);
  /// set filter for forward transform
  //int set_filters_forward(char* filtername, uint len, T* filter1, T* filter2,
  //    T* filter3 = NULL, T* filter4 = NULL);
  /// set filters for inverse transform
  //int set_filters_inverse(T* filter1, T* filter2, T* filter3 = NULL,
  //    T* filter4 = NULL);
  /// Add wavelet
  //int add_wavelet(Wavelets W, T alpha=1.0f);

 public:
  /// Image (input or result of reconstruction), on device
  T* m_image;
  /// Wavelet coefficients, on device
  std::unique_ptr<CoeffContainerT> m_coeffs;
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
