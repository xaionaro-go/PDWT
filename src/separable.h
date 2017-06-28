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

template <class Filt=void, class... Filtn>
struct MaxTapSizeLeft {
  static const constexpr int value =
    std::max(Filt::TapSizeLeft,MaxTapSizeLeft<Filtn...>::value);
};
template <> 
struct MaxTapSizeLeft<void> {
  static const constexpr int value = 0;
};
template <class Filt=void, class... Filtn>
struct MaxTapSizeRight {
  static const constexpr int value =
    std::max(Filt::TapSizeRight,MaxTapSizeRight<Filtn...>::value);
};
template <> 
struct MaxTapSizeRight<void> {
  static const constexpr int value = 0;
};
template <class Filt=void, class... Filtn>
struct MaxTapHalfSizeLeft {
  static const constexpr int value =
    std::max(Filt::TapHalfSizeLeft,MaxTapHalfSizeLeft<Filtn...>::value);
};
template <> 
struct MaxTapHalfSizeLeft<void> {
  static const constexpr int value = 0;
};
template <class Filt=void, class... Filtn>
struct MaxTapHalfSizeRight {
  static const constexpr int value =
    std::max(Filt::TapHalfSizeRight,MaxTapHalfSizeRight<Filtn...>::value);
};
template <>
struct MaxTapHalfSizeRight<void> {
  static const constexpr int value = 0;
};
template <class Filt=void, class... Filtn>
struct MaxTapHalfFloorSizeLeft {
  static const constexpr int value =
    std::max(Filt::TapHalfFloorSizeLeft,
      MaxTapHalfFloorSizeLeft<Filtn...>::value);
};
template <> 
struct MaxTapHalfFloorSizeLeft<void> {
  static const constexpr int value = 0;
};
template <class Filt=void, class... Filtn>
struct MaxTapHalfCeilSizeRight {
  static const constexpr int value =
    std::max(Filt::TapHalfCeilSizeRight,
    MaxTapHalfCeilSizeRight<Filtn...>::value);
};
template <> 
struct MaxTapHalfCeilSizeRight<void> {
  static const constexpr int value = 0;
};

template<typename T, typename I, typename J, class Filt, class... Filtn>
struct EvenSubSampledAccumulator {
  template<class In, class... InN>
  static constexpr T acc(I filtIdx, J inIdx, const In in, const InN... inn) {
   return EvenSubSampledAccumulator<T,I,J,Filt>::acc(filtIdx, inIdx, in)+
     EvenSubSampledAccumulator<T,I,J,Filtn...>::acc(filtIdx, inIdx, inn...);
  }
};
template<typename T, typename I, typename J, class Filt>
struct EvenSubSampledAccumulator<T,I,J,Filt> {
 template<class In>
 static constexpr T acc(I filtIdx, J inIdx, In in) {
    if ((filtIdx>=-Filt::TapHalfSizeLeft)&&(filtIdx<=Filt::TapHalfSizeRight)) {
	  return Filt::Buff[2*(filtIdx+Filt::TapHalfSizeLeft)+
        Filt::EvenSubSampOffset]*in[inIdx];
    } else {
      return (T)0;
    }
  }
};
template<typename T, typename I, typename J>
struct EvenSubSampledAccumulator<T,I,J,void> {
 static constexpr T acc(I index) { return (T)0; }
};

template<typename T, typename I, typename J, class Filt, class... Filtn>
struct OddSubSampledAccumulator {
 template<class In, class... InN>
 static constexpr T acc (I filtIdx, J inIdx, In in, InN... inn) {
    return  OddSubSampledAccumulator<T,I,J,Filt>::acc(filtIdx, inIdx, in) +
      OddSubSampledAccumulator<T,I,J,Filtn...>::acc(filtIdx, inIdx, inn...);
   }
};

template<typename T, typename I, typename J, class Filt>
struct OddSubSampledAccumulator<T,I,J,Filt> {
 template<class In>
 static constexpr T acc(I filtIdx, J inIdx, In in) {
    if ((filtIdx>=-Filt::TapHalfFloorSizeLeft)&&
        (filtIdx<=Filt::TapHalfCeilSizeRight)) {
	  return Filt::Buff[2*(filtIdx+Filt::TapHalfFloorSizeLeft)+
        Filt::OddSubSampOffset]*in[inIdx];
    } else {
      return (T)0;
    }
  }
};
template<typename T, typename I, typename J>
struct OddSubSampledAccumulator<T,I,J,void> {
  static constexpr T acc(I index) { return (T)0; }
};

template<typename T, typename U>
class Accumulator {
 public:
  Accumulator(T* ptr): acc(0), dst(ptr) {}
  void accumulate(T val) {
    acc+=val;
  }
  void write(size_t idx) {
    dst[idx]=acc;
  }
  void reset() {
    acc=0;
  }
protected:
  U acc;
  T* dst;
};

template<class Filt, class... Filtn>
struct Updater {
  template<typename I, typename M, class Acc, class... AccN>
  static void accumulate(I idx, M mult, Acc&& acc, AccN&& ... accn) {
    Updater<Filt>::accumulate(idx, mult, std::forward<Acc>(acc));
    Updater<Filtn...>::accumulate(idx, mult, std::forward<AccN>(accn)...);
  }
  template<typename I, class Acc, class... AccN>
  static void write(I idx, Acc&& acc, AccN&& ... accn) {
    Updater<Filt>::write(idx, std::forward<Acc>(acc));
    Updater<Filtn...>::write(idx, std::forward<AccN>(accn)...);
  }

};

/// Actually performs the accumulation
template<class Filt>
struct Updater<Filt> {
  template<typename I, typename M, class Acc>
  static void accumulate(I idx, M mult, Acc&& acc) {
    if((idx>=-Filt::TapSizeLeft) && (idx<=Filt::TapSizeRight)) {
      acc.accumulate(mult*Filt::Buff[Filt::TapSizeLeft+idx]);
    }
  }
  template<typename I, class Acc>
  static void write(I idx, Acc&& acc) {
   acc.write(idx);
  }
};

template<class Acc=void, class... AccN>
struct Resetter {
  static void reset(Acc&& acc, AccN&& ... accn) {
    acc.reset();
    Resetter<AccN...>::reset(std::forward<AccN>(accn)...);
  }
};
template<>
struct Resetter<void> {
  static void reset() {}
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
template<typename T, class... Filtn>
class SeparableSubsampledConvolutionEngine {
 public:
  /// Defaulted constructor
  SeparableSubsampledConvolutionEngine()=default;
  /// Default destructor
  virtual ~SeparableSubsampledConvolutionEngine()=default;

  /// The main method : perform Subsampled convolution on one row
  template<class... AccN>
  static int PerformSubsampledFilteringXRef(const T* in, int Nx,
      AccN&&... accn) {

    int Nx_is_odd = (Nx & 1);
    int NxOut = (Nx + Nx_is_odd)/2;

    // Loop over output image
    for (int ox=0; ox<NxOut; ox++) {
      Resetter<AccN...>::reset(std::forward<AccN>(accn)...);
      // Loop over filter size, with periodic boundary conditions
      // TODO TN: this loop can actually be written as a compile time loop
      // #pragma unroll Filt::TapSize
      for (int fx=-MaxTapSizeLeft<Filtn...>::value;
          fx<=MaxTapSizeRight<Filtn...>::value; fx++) {
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
        // Update each buffer with its respective filter
        Updater<Filtn...>::accumulate(fx, in[ix], std::forward<AccN>(accn)...);
      }
      Updater<Filtn...>::write(ox,std::forward<AccN>(accn)...);
    }
    return 1;
  }
};

/** \class SeparableUpsampledConvolutionEngine
 * \brief Code for the separable upsampled convolution.
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class... Filtn>
class SeparableUpsampledConvolutionEngine {
 public:
  /// Defaulted constructor
  SeparableUpsampledConvolutionEngine()=default;
  /// Default destructor
  virtual ~SeparableUpsampledConvolutionEngine()=default;

  /// The main method : perform Subsampled convolution on one row
  template<class... InN>
  static int PerformUpsampledFilteringXRef(int NxIn, int NxOut, T* out,
    InN... inn) {

    // Loop over output image
    for (int lox=0; lox<NxOut; lox++) {
      int max_x = NxIn-1;
      int ixCentral = lox/2;
      T acc = (T)0;

      if ((lox%2)==0) {
	    //TODO TN Filter loop, can be turned into an explicit compile time loop
		for (int fx=-MaxTapHalfSizeLeft<Filtn...>::value;
            fx<=MaxTapHalfSizeRight<Filtn...>::value;fx++) {
		  int idx_x = ixCentral + fx;
		  if (idx_x<0) {
			idx_x += NxIn;
		  }
		  if (idx_x>max_x) {
			idx_x -= NxIn;
		  }
          acc+=EvenSubSampledAccumulator<T,int,int,Filtn...>::acc(
            fx, idx_x, inn...);
          /*if ((fx>=-FiltLow::TapHalfSizeLeft)&&
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
		  }*/
		}
      } else {
	    //TODO TN Filter loop, can be turned into an explicit compile time loop
		for (int fx=-MaxTapHalfFloorSizeLeft<Filtn...>::value;
            fx<=MaxTapHalfCeilSizeRight<Filtn...>::value; fx++) {
		  int idx_x = ixCentral + fx;
		  if (idx_x<0) {
			idx_x += NxIn;
		  }
		  if (idx_x>max_x) {
			idx_x -= NxIn;
		  }
	      acc+=OddSubSampledAccumulator<T,int,int,Filtn...>::acc(
            fx, idx_x, inn...);
	      // conditional update with respective filter
		  /*if ((fx>=-FiltLow::TapHalfFloorSizeLeft)&&
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
		  }*/
		}

      }
      // Update each buffer with its respective filter
      out[lox] = acc;
    }
    return 1;
  }
};



#endif //SEPARABLE_H
