// STL
#include <cassert>

// Pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

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
      throw std::runtime_error("Number of dimensions must be one");
    }
    T* ptr = (double *) buffer.ptr;
    size_t size = buffer.size;

    if (m_name=="Daub2") {
      m_pWavelet = std::make_unique<Daub2_1D<T>>(
        ptr,size,1,1,m_doCycleSpinning,m_name,m_nbLevel);
    } else if (m_name=="dtwAnto97QSHIFT6") {
      m_pWavelet = std::make_unique<dtwAnto97QSHIFT6_1D<T>>(
        ptr,size,1,1,m_doCycleSpinning,m_name,m_nbLevel);
      //
    } else {
      assert(false);
    }
  }
 protected:
  int m_nbLevel;
  bool m_doCycleSpinning;
  std::string m_name;
  std::unique_ptr<WaveletWrapper> m_pWavelet;
};

PYBIND11_PLUGIN(Wavelet1D) {
  py::module m("Wavelet1D", "pybind11 wavelet binding");
  py::class_<Wavelet1DWrapper<float>>(m, "Wavelet1D",py::dynamic_attr())
  .def(py::init<int,bool,const std::string &>())
  .def("Initialize", &Wavelet1DWrapper<float>::Initialize);
  return m.ptr();
}

