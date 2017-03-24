#ifndef WAVELET_H
#define WAVELET_H

// STL
#include <list>
#include <memory>
#include <string>

// Local

/** \struct w_info
 * \brief Description of the workload
 *
 * \author Pierre Paleo
 */
struct w_info {
  /// Number of columns in the image
  int Nc;
  /// Number of rows of the image (2D)
  int Nr;
  /// Number of slices in the image (3D)
  int Ns;
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
  /// Error when creating the Wavelet instance
  W_CREATION_ERROR,
  /// Error when computing the forward transform
  W_FORWARD_ERROR,
  /// Error when computing the inverse transform
  W_INVERSE_ERROR,
  /// Error when thresholding the coefficients
  W_THRESHOLD_ERROR
};

/** \class Wavelet
 * \brief The wavelet class manage the lifecycle of the wavelet transform
 * coefficients.
 *
 * \author Pierre Paleo
 */
template<typename T, class CoeffContainerT, class WaveletSchemeT>
class Wavelet {
 public:
  /// Constructor with zero initialization
  Wavelet();
  /// Constructor : Wavelet from image
  Wavelet(T* img, int Nc, int Nr, int Ns, bool doCycleSpinning,
    const std::string& wname, int level); 
  /// Copy constructor
  Wavelet(const Wavelet& w)=delete;
  /// Default destructor
  virtual ~Wavelet()=default;
  
  /// Print wavelet transform informations
  virtual void print_informations();

  /// Forward wavelet tranform
  virtual int forward();
  /// Backward wavelet transform: transpose of the forward transpose
  virtual int backward();
  /// Inverse of the wavelet tranform
  virtual int inverse();
 
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
  virtual T* get_image(T* img);
    /// Return a pointer to wavelet coefficients
  virtual T* get_coeff();
  /// Set input image to be transformed in the wavelet domain
  int virtual set_image(T* img);
  /// Set coefficients to be reconstructed into an image
  int virtual set_coeff(T* coeff);
  /// set filter for forward transform
  //int set_filters_forward(char* filtername, uint len, T* filter1, T* filter2,
  //    T* filter3 = NULL, T* filter4 = NULL);
  /// set filters for inverse transform
  //int set_filters_inverse(T* filter1, T* filter2, T* filter3 = NULL,
  //    T* filter4 = NULL);
  /// Add wavelet
  //int add_wavelet(Wavelet W, T alpha=1.0f);

 public:
  /// Image (input or result of reconstruction), on device
  T* m_image;
  /// Wavelet coefficients, on device
  std::unique_ptr<CoeffContainerT> m_coeff;
  /** If cycle spinning is enabled, use image shifting for translation
   * invariant denoising
   */
  bool m_doCycleSpinning;
  /// Current shift for the cycle spinning process
  std::list<int> m_currentShift;
  /// Wavelet name
  std::string m_name;
  /// Informations about the wavelet tranform
  w_info m_info;
  /// Current state of the wavelet tranform
  w_state m_state;
  /// Number of level od dyadic decomposition
  int m_level;
};


#endif //WAVELET_H