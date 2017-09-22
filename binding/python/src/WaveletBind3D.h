// STL
#include <cassert>

// Local
#include "filters.h"
#include "Wavelet3D.h"

namespace py = pybind11;

/**
 * One may be interested in taking a look at Wavelet3D.cpp in the src
 * directory in order to instanciate the chosen wavelet transform with
 * an existing c++ version
 */

template<typename T>
class Wavelet3DWrapper {
 public:
  Wavelet3DWrapper(int nbLevel=1, bool doCycleSpinning=false,
      const std::string &name="Daub4") : m_nbLevel(nbLevel), 
      m_doCycleSpinning(doCycleSpinning), m_name(name) {}

  void Initialize(py::array_t<T> image) {
    auto buffer = image.request();
    if (buffer.ndim != 3) {
      throw std::runtime_error("Wavelet3DWrapper::Initialize : "
        "Number of dimensions must be one");
    }
    T* ptr = static_cast<T *>(buffer.ptr);
    size_t size = buffer.size;

    if (m_name=="Daub2") {
      m_pWavelet = std::make_unique<Daub2_3D<T>>(
        ptr,buffer.shape[0],buffer.shape[1],buffer.shape[2],
        m_doCycleSpinning,m_name,m_nbLevel);
    } else if (m_name=="Daub3") {
      m_pWavelet = std::make_unique<Daub3_3D<T>>(
        ptr,buffer.shape[0],buffer.shape[1],buffer.shape[2],
        m_doCycleSpinning,m_name,m_nbLevel);
    } else if (m_name=="Daub4") {
      m_pWavelet = std::make_unique<Daub4_3D<T>>(
        ptr,buffer.shape[0],buffer.shape[1],buffer.shape[2],
        m_doCycleSpinning,m_name,m_nbLevel);
    } else if (m_name=="Daub5") {
      m_pWavelet = std::make_unique<Daub5_3D<T>>(
        ptr,buffer.shape[0],buffer.shape[1],buffer.shape[2],
        m_doCycleSpinning,m_name,m_nbLevel);
    } else if (m_name=="dtwAnto97QSHIFT6") {
       m_pWavelet = std::make_unique<dtwAnto97QSHIFT6_3D<T>>(
        ptr, buffer.shape[0],buffer.shape[1],buffer.shape[2],
        m_doCycleSpinning,m_name,m_nbLevel);
    } else {
      throw std::runtime_error("Wavelet3DWrapper::Initialize : "
        "Unsupported wavelet type");
    }
  }
  void forward() {
    if (m_pWavelet->forward()<0) {
      throw std::runtime_error("Wavelet3DWrapper::forward : "
        "Runtime error");
    }
  }
  void backward() {
    if (m_pWavelet->backward()<0) {
      throw std::runtime_error("Wavelet3DWrapper::backward : "
        "Runtime error");
    }
  }
  void inverse() {
    if (m_pWavelet->inverse()<0) {
      throw std::runtime_error("Wavelet3DWrapper::inverse : "
        "Runtime error");
    }
  }
  void set_image(py::array_t<T> image) {
    auto buffer = image.request();
    if (buffer.ndim != 2) {
      throw std::runtime_error("Wavelet3DWrapper::set_image : "
        "Number of dimensions must be one");
    }
    T* ptr = static_cast<T *>(buffer.ptr);
    if (m_pWavelet->set_image(ptr)<0) {
      throw std::runtime_error("Wavelet3DWrapper::set_image : "
        "Runtime error");
    }
  }
 protected:
  int m_nbLevel;
  bool m_doCycleSpinning;
  std::string m_name;
  std::unique_ptr<WaveletWrapper<T>> m_pWavelet;
};

