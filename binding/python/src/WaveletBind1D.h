// STL
#include <cassert>

// Local
#include "filters.h"
#include "Wavelet1D.h"

namespace py = pybind11;

/**
 * One may be interested in taking a look at Wavelet1D.cpp in the src
 * directory in order to instanciate the chosen wavelet transform with
 * an existing c++ version
 */

template<typename T>
class Wavelet1DWrapper {
 public:
  Wavelet1DWrapper(int nbLevel=1, bool doCycleSpinning=false,
      const std::string &name="Daub4") : m_nbLevel(nbLevel), 
      m_doCycleSpinning(doCycleSpinning), m_name(name) {}

  void Initialize(py::array_t<T> image) {
    auto buffer = image.request();
    if (buffer.ndim != 1) {
      throw std::runtime_error("Wavelet1DWrapper::Initialize : "
        "Number of dimensions must be one");
    }
    T* ptr = static_cast<T *>(buffer.ptr);
    size_t size = buffer.size;

    if (m_name=="Daub2") {
      m_pWavelet = std::make_unique<Daub2_1D<T>>(
        ptr,size,1,1,m_doCycleSpinning,m_name,m_nbLevel);
    } else if (m_name=="dtwAnto97QSHIFT6") {
      m_pWavelet = std::make_unique<dtwAnto97QSHIFT6_1D<T>>(
        ptr,size,1,1,m_doCycleSpinning,m_name,m_nbLevel);
    } else {
      throw std::runtime_error("Wavelet1DWrapper::Initialize : "
        "Unsupported wavelet type");
    }
  }
  void forward() {
    if (m_pWavelet->forward()<0) {
      throw std::runtime_error("Wavelet1DWrapper::forward : "
        "Runtime error");
    }
  }
  void backward() {
    if (m_pWavelet->backward()<0) {
      throw std::runtime_error("Wavelet1DWrapper::backward : "
        "Runtime error");
    }
  }
  void inverse() {
    if (m_pWavelet->inverse()<0) {
      throw std::runtime_error("Wavelet1DWrapper::inverse : "
        "Runtime error");
    }
  }
  void set_image(py::array_t<T> image) {
    auto buffer = image.request();
    if (buffer.ndim != 1) {
      throw std::runtime_error("Wavelet1DWrapper::set_image : "
        "Number of dimensions must be one");
    }
    T* ptr = static_cast<T *>(buffer.ptr);
    if (m_pWavelet->set_image(ptr)<0) {
      throw std::runtime_error("Wavelet1DWrapper::set_image : "
        "Runtime error");
    }
  }
/*  py::array_t<T> get_coeff() const {
    py::array_t<T>();
    auto buffer = image.request();
    if (buffer.ndim != 1) {
      throw std::runtime_error("Wavelet1DWrapper::get_image : "
        "Number of dimensions must be one");
    }
    T* ptr = static_cast<T*>(buffer.ptr);
    if (m_pWavelet->get_image(ptr)<0) {
      throw std::runtime_error("Wavelet1DWrapper::get_image : "
        "Runtime error");
    }
  }*/
 protected:
  int m_nbLevel;
  bool m_doCycleSpinning;
  std::string m_name;
  std::unique_ptr<WaveletWrapper<T>> m_pWavelet;
};

