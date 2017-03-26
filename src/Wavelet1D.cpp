// Local
#include "Wavelet1D.h"

// Local
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

// Instanciating all 1d wavelet types at once
template class DB1DWt<float>;
template class DB1DWt<double>;
