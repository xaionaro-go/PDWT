

// Local
#include "filters.h"
#include <Wavelet1D.h>

/**
 * One may be interested in taking a look at Wavelet1D.cpp in the src
 * directory in order to instanciate the chosen wavelet transform with
 * an existing c++ version
 */


template<> Wavelet1D<float,PackedContainer1D<float>,Daub2<float>>;
