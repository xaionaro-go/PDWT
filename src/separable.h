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
 * \brief Code for the separable upsampled convolution.
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class FiltLow, class FiltHigh>
class SeparableUpsampledConvolutionEngine {
 public:
  /// Defaulted constructor
  SeparableUpsampledConvolutionEngine()=default;
  /// Default destructor
  virtual ~SeparableUpsampledConvolutionEngine()=default;

  /// The main method : perform Subsampled convolution on one row
  static int PerformUpsampledFilteringXRef(const T* inLow, T* inHigh, int NxIn,
    int NxOut, T* out) {

    // Loop over output image
    for (int lox=0; lox<NxOut; lox++) {
      int ox = lox + FiltLow::IsHalfSizeOdd?0:1;
      int ixCentral = ox/2;
      int max_x = NxIn-1;
      //si index impair: pas d'offset, sinon offset 1
	  int offset_x = 1-(ox&1);
      T acc = (T)0;
 
      if (offset_x==0) {
		//TODO TN Filter loop, can be turned into an explicit compile time loop
		for (int jx = 0; jx <= FiltLow::TapHalfSize; jx++) {
			int idx_x = ixCentral - FiltLow::TapHalfSizeLeft + jx;
			if (idx_x<0) idx_x += NxIn;
			if (idx_x>max_x) idx_x -= NxIn;
			//int fAddr = FiltLow::TapSize -1 - (2*jx + offset_x);
			int fAddr = FiltLow::TapSize -1 - (2*jx);
			// Update each buffer with its respective filter
			acc += inLow[idx_x] * FiltLow::Buff[fAddr];
			acc += inHigh[idx_x] * FiltHigh::Buff[fAddr];
			//Updater<FiltLow>::update(fAddr, inLow[idx_x], out+ox);
		}
      } else {
		//TODO TN Filter loop, can be turned into an explicit compile time loop
		for (int jx = 0; jx <= FiltLow::TapHalfSize; jx++) {
			int idx_x = ixCentral - FiltLow::TapHalfSizeLeft + jx;
			if (idx_x<0) idx_x += NxIn;
			if (idx_x>max_x) idx_x -= NxIn;
			//int fAddr = FiltLow::TapSize -1 - (2*jx + offset_x);
			int fAddr = FiltLow::TapSize -1 - (2*jx+1);
			// Update each buffer with its respective filter
			acc += inLow[idx_x] * FiltLow::Buff[fAddr];
			acc += inHigh[idx_x] * FiltHigh::Buff[fAddr];
			//Updater<FiltLow>::update(fAddr, inLow[idx_x], out+ox);
		}
      }
      // Update each buffer with its respective filter
      out[lox] = acc;
    }
    return 1;
  }
};


// must be run with grid size = (2*Nr, 2*Nc) ; Nc = numcols of input coeffs. Here the param Nr is actually doubled wrt Nr_coeffs because of the vertical oversampling.
// pass 2 : (tmp1, tmp2)  ==> Horiz convol with IL, IH  + horiz oversampling ==> I
/**
One has to keep in mind that forward transform in a convolution THEN subsampling
ex:
|1 0 0 0| |1 -1 0 -1|
|0 0 1 0| |-1 1 -1 0|
          |0 -1 1 -1|
          |-1 0 -1 1|

Transpose of that is then:

|1 -1 0 -1| |1 0|
|-1 1 -1 0| |0 0|
|0 -1 1 -1| |0 1|
|-1 0 -1 1| |0 0|
So that only one out of 2 coefficients in the convolution is useful

__global__ void w_kern_inverse_pass2(
	DTYPE* tmp1, //input lowfrequency space
    DTYPE* tmp2, //input highfrequency space
	DTYPE* img,  //output image: upsampled in the x direction
    int Nr,      //in/out size in y
	int Nc,      //input size in x
	int Nc2,     //output size in x ater upsampling
	int hlen     //size of the filter
  ) {

    int gidx = threadIdx.x + blockIdx.x*blockDim.x;  //gidx: output address x (up to twice the size of the input)
    int gidy = threadIdx.y + blockIdx.y*blockDim.y;  //gidy: input/output  address y

    if (gidy < Nr && gidx < Nc2) { // horiz oversampling : Input (Nr*2, Nc) => Output (Nr*2, Nc*2)

        int hL, hR;
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
        int max_x = Nc-1;
        int offset_x = 1-(gidx & 1);i//si index impair: pas d'offset, sinon offset 1

        for (int jx = 0; jx <= hR+hL; jx++) {
            int idx_x = jx1 + jx;
            if (idx_x<0) idx_x += Nc;
            if (idx_x>max_x) idx_x -= Nc;
            res_1 += tmp1[gidy*Nc + idx_x] * c_kern_IL[hlen-1 - (2*jx + offset_x)];
            res_2 += tmp2[gidy*Nc + idx_x] * c_kern_IH[hlen-1 - (2*jx + offset_x)];
        }
        //Si half kernel est impair: on peut ecrire avec le bon mapping
        if ((hlen2 & 1) == 1) img[gidy * Nc2 + gidx] = res_1 + res_2;
        //Sinon on decale l'ecriture de 1 sur la gauche
        else img[gidy * Nc2 + (gidx-1)] = res_1 + res_2;
    }
}
int jx1 = c - gidx/2;//= -(gidx/2-hl)
        int jx2 = Nc - 1 - gidx/2 + c;//=Nc-1:adresse dernier elem entree, Nc-1-(gidx/2-hl)
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
*/
#endif //SEPARABLE_H
