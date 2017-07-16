// STL
#include <cassert>

// Pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// Local
#include "WaveletBind1D.h"
#include "WaveletBind2D.h"

namespace py = pybind11;

PYBIND11_MODULE(pyPDWT, m) {
  m.doc() = "pyPDWT : pybind11 wavelet binding";

  py::class_<Wavelet1DWrapper<float>> Wavelet1D(m, "Wavelet1D",
    py::dynamic_attr());
  Wavelet1D
    .def(py::init<int,bool,const std::string &>())
    .def("Initialize", &Wavelet1DWrapper<float>::Initialize)
    .def("forward", &Wavelet1DWrapper<float>::forward)
    .def("backward", &Wavelet1DWrapper<float>::backward)
    .def("inverse", &Wavelet1DWrapper<float>::inverse)
    .def("set_image", &Wavelet1DWrapper<float>::set_image);

  py::class_<Wavelet2DWrapper<float>> Wavelet2D(m, "Wavelet2D",
    py::dynamic_attr());
  Wavelet2D
    .def(py::init<int,bool,const std::string &>())
    .def("Initialize", &Wavelet2DWrapper<float>::Initialize)
    .def("forward", &Wavelet2DWrapper<float>::forward)
    .def("backward", &Wavelet2DWrapper<float>::backward)
    .def("inverse", &Wavelet2DWrapper<float>::inverse)
    .def("set_image", &Wavelet2DWrapper<float>::set_image);
}

