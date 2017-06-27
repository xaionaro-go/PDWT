#ifndef SEPARABLE_H
#define SEPARABLE_H

// STL

// Local
#include "filters.h"

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
*/

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
//TODO TN: replace second line by first one
//template<typename T, class Filt, class... Filtn>
template<typename T, class FiltLow, class FiltHigh>
class SeparableSubsampledConvolutionEngine {
 public:
  /// Defaulted constructor
  SeparableSubsampledConvolutionEngine()=default;
  /// Default destructor
  virtual ~SeparableSubsampledConvolutionEngine()=default;

  /// The main method : perform Subsampled convolution on one row
  //template<typename... O>
  //TODO TN: to be modified
  static int PerformSubsampledFilteringXRef(const T* in, int Nx,
      T* outlow, T* outhigh) {
//      O... out) {

    int Nx_is_odd = (Nx & 1);
    int NxOut = (Nx + Nx_is_odd)/2;

    // Loop over output image
    for (int ox=0; ox<NxOut; ox++) {
      //TOD TN: delete next 2 lines
      outlow[ox]=0;
      outhigh[ox]=0;
      // Loop over filter size, with periodic boundary conditions
      // TODO TN: this loop can actually be written as a compile time loop
      // #pragma unroll Filt::TapSize
      for (int fx=-std::max(FiltLow::TapSizeLeft,FiltHigh::TapSizeLeft);
          fx<=std::max(FiltLow::TapSizeRight,FiltHigh::TapSizeRight); fx++) {
        int ix = ox*2 + fx;
        // if N is odd, image is virtually extended
        if (ix < 0) {
          ix += (Nx + Nx_is_odd);
        }
        // no "else if", since idx_x can be > N-1  after being incremented
        if (ix > Nx-1) {
          // if N is odd, repeat the right-most element
          if ((ix == Nx) && (Nx_is_odd)) {
            ix--;
          } else {
           // if N is odd, image is virtually extended
           ix -= (Nx + Nx_is_odd);
          }
        }
        //std::cout<<"Accumulate low image idx "<<ix<<": "<<in[ix]
        //  <<" x "<<FiltLow::Buff[fx+FiltLow::TapSizeLeft]
        //  <<" filt idx: "<<fx+FiltLow::TapSizeLeft<<std::endl;

        // Update each buffer with its respective filter
        //Updater<Filt,Filtn...>::update(fx, in[ix], out+ox...);
        
        //actualy, this is a conditional update
        if ((fx>=-FiltLow::TapSizeLeft)&&(fx<=FiltLow::TapSizeRight)) {
          outlow[ox]+=in[ix]*FiltLow::Buff[fx+FiltLow::TapSizeLeft];
        }
        if ((fx>=-FiltHigh::TapSizeLeft)&&(fx<=FiltHigh::TapSizeRight)) {
          outhigh[ox]+=in[ix]*FiltHigh::Buff[fx+FiltHigh::TapSizeLeft];
        }
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
      int max_x = NxIn-1;
      int ixCentral = lox/2;
      T acc = (T)0;

      if ((lox%2)==0) {
	    //TODO TN Filter loop, can be turned into an explicit compile time loop
		for (int 
			fx=-std::max(FiltLow::TapHalfSizeLeft,FiltHigh::TapHalfSizeLeft);
			fx<=std::max(FiltLow::TapHalfSizeRight,FiltHigh::TapHalfSizeRight);
			fx++) {
		  int idx_x = ixCentral + fx;
		  if (idx_x<0) {
			idx_x += NxIn;
		  }
		  if (idx_x>max_x) {
			idx_x -= NxIn;
		  }
		  // conditional update with respective filter
		  if ((fx>=-FiltLow::TapHalfSizeLeft)&&
			 (fx<=FiltLow::TapHalfSizeRight)) {
		   int fAddr = 2*(fx+FiltLow::TapHalfSizeLeft)+
             FiltLow::EvenSubSampOffset;
		   acc += inLow[idx_x] * FiltLow::Buff[fAddr];
		  }
		  if ((fx>=-FiltHigh::TapHalfSizeLeft)&&
			  (fx<=FiltHigh::TapHalfSizeRight)) {
			int fAddr = 2*(fx+FiltHigh::TapHalfSizeLeft)+
              FiltHigh::EvenSubSampOffset;
			acc += inHigh[idx_x] * FiltHigh::Buff[fAddr];
		  }
		}
      } else {
	    //TODO TN Filter loop, can be turned into an explicit compile time loop
		for (int 
			fx=-std::max(FiltLow::TapHalfFloorSizeLeft,
              FiltHigh::TapHalfFloorSizeLeft);
			fx<=std::max(FiltLow::TapHalfCeilSizeRight,
              FiltHigh::TapHalfCeilSizeRight);
			fx++) {
		  int idx_x = ixCentral + fx;
		  if (idx_x<0) {
			idx_x += NxIn;
		  }
		  if (idx_x>max_x) {
			idx_x -= NxIn;
		  }
		  // conditional update with respective filter
		  if ((fx>=-FiltLow::TapHalfFloorSizeLeft)&&
		      (fx<=FiltLow::TapHalfCeilSizeRight)) {
	   	    int fAddr = 2*(fx+FiltLow::TapHalfFloorSizeLeft)+
              FiltLow::OddSubSampOffset;
            acc += inLow[idx_x] * FiltLow::Buff[fAddr];
		  }
		  if ((fx>=-FiltHigh::TapHalfFloorSizeLeft)&&
			  (fx<=FiltHigh::TapHalfCeilSizeRight)) {
	        int fAddr = 2*(fx+FiltHigh::TapHalfFloorSizeLeft)+
              FiltHigh::OddSubSampOffset;
		    acc += inHigh[idx_x] * FiltHigh::Buff[fAddr];
		  }
		}

      }
      // Update each buffer with its respective filter
      out[lox] = acc;
    }
    return 1;
  }
};



#endif //SEPARABLE_H
