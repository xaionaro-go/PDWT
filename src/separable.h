#ifndef SEPARABLE_H
#define SEPARABLE_H

// STL

// Local
#include "filters.h"


// TODO TN: you should replace this updater by a slightly more complex struct
// which is also in charge of storing an accumulator, and has actually two
// methods: accumulate and update
template<class Filt, class... Filtn>
struct Updater {
  template<typename I, typename M, typename O, typename... On>
  static void update(I idx, M mult, O out, On... outn) {
    Updater<Filt>::update(idx, mult, out);
    Updater<Filtn...>::update(idx, mult, outn...);
  }
};

/// Actually performs the accumulation
template<class Filt>
struct Updater<Filt> {
  template<typename I, typename M, typename O>
  static void update(I idx, M mult, O out) {
    *out+=mult*Filt::Buff[idx+Filt::TapSizeLeft];
  }
};

/** \class SeparableSubsampledConvolutionEngine
 * \brief Code for the separable subsampled convolution. This class is a
 * variadic template class, because it can handle multiple filtering for each
 * main loop, assuming the filters have the same size
 *
 * TODO TN: perf issue: you should use temporary for accumulation
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class Filt, class... Filtn>
class SeparableSubsampledConvolutionEngine {
 public:
  /// Defaulted constructor
  SeparableSubsampledConvolutionEngine()=default;
  /// Default destructor
  virtual ~SeparableSubsampledConvolutionEngine()=default;

  /// The main method : perform Subsampled convolution on one row
  template<typename... O>
  static int PerformSubsampledFilteringXRef(const T* in, int Nx,
      O... out) {

    int Nx_is_odd = (Nx & 1);
    int NxOut = (Nx + Nx_is_odd)/2;

    // Loop over output image
    for (int ox=0; ox<NxOut; ox++) {
      // Loop over filter size, with periodic boundary conditions
      // TODO TN: this loop can actually be written as a compile time loop
      // #pragma unroll Filt::TapSize
      for (int fx=-Filt::TapSizeLeft; fx<Filt::TapSizeRight; fx++) {
        int ix = ox*2 + fx;
        // if N is odd, image is virtually extended
        if (ix < 0) ix += (Nx + Nx_is_odd);
        // no "else if", since idx_x can be > N-1  after being incremented
        if (ix > Nx-1) {
          // if N is odd, repeat the right-most element
          if ((ix == Nx) && (Nx_is_odd))
            ix--;
          // if N is odd, image is virtually extended
          else
            ix -= (Nx + Nx_is_odd);
        }
        // Update each buffer with its respective filter
        Updater<Filt,Filtn...>::update(fx, in[ix], out+ox...);
      }
    }
    return 1;
  }
};

/** \class SeparableUpsampledConvolutionEngine
 * \brief Code for the separable upsampled convolution. This class is a
 * variadic template class, because it can handle multiple filtering for each
 * main loop, assuming the filters have the same size
 *
 * TODO TN: perf issue: you should use temporary for accumulation
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class Filt, class... Filtn>
class SeparableUpsampledConvolutionEngine {
 public:
  /// Defaulted constructor
  SeparableUpsampledConvolutionEngine()=default;
  /// Default destructor
  virtual ~SeparableUpsampledConvolutionEngine()=default;

  /// The main method : perform Subsampled convolution on one row
  template<typename... O>
  static int PerformUpsampledFilteringXRef(const T* in, int Nx,
      O... out) {

    int Nx_is_odd = (Nx & 1);
    int NxOut = (Nx + Nx_is_odd)/2;

    // Loop over output image
    for (int ox=0; ox<NxOut; ox++) {
      // Loop over filter size, with periodic boundary conditions
      // TODO TN: this loop can actually be written as a compile time loop
      // #pragma unroll Filt::TapSize
      for (int fx=-Filt::TapSizeLeft; fx<Filt::TapSizeRight; fx++) {
        int ix = ox*2 + fx;
        // if N is odd, image is virtually extended
        if (ix < 0) ix += (Nx + Nx_is_odd);
        // no "else if", since idx_x can be > N-1  after being incremented
        if (ix > Nx-1) {
          // if N is odd, repeat the right-most element
          if ((ix == Nx) && (Nx_is_odd))
            ix--;
          // if N is odd, image is virtually extended
          else
            ix -= (Nx + Nx_is_odd);
        }
        // Update each buffer with its respective filter
        Updater<Filt,Filtn...>::update(fx, in[ix], out+ox...);
      }
    }
    return 1;
  }
};

// must be run with grid size = (2*Nr, 2*Nc) ; Nc = numcols of input coeffs. Here the param Nr is actually doubled wrt Nr_coeffs because of the vertical oversampling.
// pass 2 : (tmp1, tmp2)  ==> Horiz convol with IL, IH  + horiz oversampling ==> I
__global__ void w_kern_inverse_pass2(DTYPE* tmp1, DTYPE* tmp2, DTYPE* img, int Nr, int Nc, int Nc2, int hlen) {
    int gidx = threadIdx.x + blockIdx.x*blockDim.x;
    int gidy = threadIdx.y + blockIdx.y*blockDim.y;
    if (gidy < Nr && gidx < Nc2) { // horiz oversampling : Input (Nr*2, Nc) => Output (Nr*2, Nc*2)
        int c, hL, hR;
        int hlen2 = hlen/2; // Convolutions with even/odd indices of the kernels
        if (hlen2 & 1) { // odd half-kernel size
            c = hlen2/2;
            hL = c;
            hR = c;
        }
        else { // even half-kernel size : center is shifted to the RIGHT for reconstruction.
            c = hlen2/2 - 0;
            hL = c;
            hR = c-1;
            // virtual id for shift
            // TODO : for the very first convolution (on the edges), this is not exactly accurate (?)
            gidx += 1;
        }
        //Ce passage est important: il y a autant de thread que de pixels de sortie
        //gidx/2 represente l'indice de l'image d'entree qui va devoir etre considere
        //gidx/2-c represente l'indice a partir duquel la convolution devrait commencer
        //
        int jx1 = c - gidx/2;
        int jx2 = Nc - 1 - gidx/2 + c;
        int offset_x = 1-(gidx & 1);

        DTYPE res_1 = 0, res_2 = 0;
        for (int jx = 0; jx <= hR+hL; jx++) {
            int idx_x = gidx/2 - c + jx;
            if (jx < jx1) idx_x += Nc;
            if (jx > jx2) idx_x -= Nc;

            res_1 += tmp1[gidy*Nc + idx_x] * c_kern_IL[hlen-1 - (2*jx + offset_x)];
            res_2 += tmp2[gidy*Nc + idx_x] * c_kern_IH[hlen-1 - (2*jx + offset_x)];
        }
        if ((hlen2 & 1) == 1) img[gidy * Nc2 + gidx] = res_1 + res_2;
        else img[gidy * Nc2 + (gidx-1)] = res_1 + res_2;
    }
}

#endif //SEPARABLE_H
