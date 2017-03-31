

// Local
#include "CoeffContainer1D.h"

py::class_<Matrix>(m, "Coeff1D", py::buffer_protocol())
.def_buffer( [](Matrix &m) -> py::buffer_info {
    return py::buffer_info(
      m.data(),       // Pointer to buffer
      sizeof(float),  // Size of one scalar
      py::format_descriptor<float>::format(), //Python struct-style format desc
      1,              //Number of dimension
      { m.rows(), m.cols() },  // Buffer dimension in number of elements
      { sizeof(float) * m.rows(), sizeof(float) } //Stride in byte for each idx
    );
  }
);
