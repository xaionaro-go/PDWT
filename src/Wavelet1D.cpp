// Local
#include "Wavelet1D.h"

// Local
#include "filters.h"
#include "separable.h"

template<typename T, class CoeffContainerT, class WaveletSchemeT>
Wavelet1D<T,CoeffContainerT,WaveletSchemeT>::Wavelet1D() {

}

template<typename T, class CoeffContainerT, class WaveletSchemeT>
Wavelet1D<T,CoeffContainerT,WaveletSchemeT>:: Wavelet1D(
    T* img, int Nc, int Nr, int Ns, bool doCycleSpinning,
    const std::string& wname, int level) :
    Wavelet<T,CoeffContainerT,WaveletSchemeT>(img, Nc, Nr, Ns, doCycleSpinning,
    wname, level) {

}

template<typename T, class CoeffContainerT, class WaveletSchemeT>
int Wavelet1D<T,CoeffContainerT,WaveletSchemeT>::forward() {
  return 1;
}

template<typename T, class CoeffContainerT, class WaveletSchemeT>
int Wavelet1D<T,CoeffContainerT,WaveletSchemeT>::backward() {
  return 1;
}

template<typename T, class CoeffContainerT, class WaveletSchemeT>
int Wavelet1D<T,CoeffContainerT,WaveletSchemeT>::inverse() {
  return 1;
}

// Instanciation
template class Wavelet1D<float,PackedContainer1D<float>,Daub2<float>>;
template class Wavelet1D<double,PackedContainer1D<double>,Daub2<double>>;
