#ifndef SEPARABLE_H
#define SEPARABLE_H

// STL
#include<tuple>

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
struct EvenSubsampledAccumulator {
  template<class In, class... InN>
  static constexpr T acc(I filtIdx, J inIdx, const In in, const InN... inn) {
   return EvenSubsampledAccumulator<T,I,J,Filt>::acc(filtIdx, inIdx, in)+
     EvenSubsampledAccumulator<T,I,J,Filtn...>::acc(filtIdx, inIdx, inn...);
  }
};
template<typename T, typename I, typename J, class Filt>
struct EvenSubsampledAccumulator<T,I,J,Filt> {
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
struct EvenSubsampledAccumulator<T,I,J,void> {
 static constexpr T acc(I index) { return (T)0; }
};

template<typename T, typename I, typename J, class Filt, class... Filtn>
struct OddSubsampledAccumulator {
 template<class In, class... InN>
 static constexpr T acc (I filtIdx, J inIdx, In in, InN... inn) {
    return  OddSubsampledAccumulator<T,I,J,Filt>::acc(filtIdx, inIdx, in) +
      OddSubsampledAccumulator<T,I,J,Filtn...>::acc(filtIdx, inIdx, inn...);
   }
};

template<typename T, typename I, typename J, class Filt>
struct OddSubsampledAccumulator<T,I,J,Filt> {
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
struct OddSubsampledAccumulator<T,I,J,void> {
  static constexpr T acc(I index) { return (T)0; }
};

template<typename T, typename U, typename V, typename I>
class Accumulator {
 public:
  Accumulator(T* srcPtr, V* dstPtr, I srcStride=1, I dstStride=1):
    m_acc(0), m_srcPtr(srcPtr), m_dstPtr(dstPtr),
    m_srcStride(srcStride), m_dstStride(dstStride)  {}

  void accumulate(I srcIdx, U filt) {
    m_acc+=m_srcPtr[srcIdx*m_srcStride]*filt;
  }
  
  void write(I dstIdx) {
    m_dstPtr[dstIdx*m_dstStride]=m_acc;
  }
  void reset() {
    m_acc=0;
  }
  void incrementSrcDstPtr(I srcInc, I dstInc) {
    m_srcPtr+=srcInc;
    m_dstPtr+=dstInc;
  }
  void incrementSrcDstPtrStrided(I srcInc, I dstInc) {
    m_srcPtr+=srcInc*m_srcStride;
    m_dstPtr+=dstInc*m_dstStride;
  }
 protected:
  T* m_srcPtr;
  U m_acc;
  V* m_dstPtr;
  I m_srcStride;
  I m_dstStride;
};

template<typename T, typename U, typename I, typename J, class... Filtn>
class SubsampledAccumulator {
 public:
  SubsampledAccumulator(T* dstPtr, I srcStride=1, I dstStride=1):
    m_dstPtr(dstPtr), m_acc(0), m_srcStride(srcStride),
    m_dstStride(dstStride)  {}

  template<class... InN>
  void EvenSubsampledAccumulate(I filtIdx, J inIdx, const InN... inn) {
    m_acc+=EvenSubsampledAccumulator<T,I,J,Filtn...>::acc(
      filtIdx, inIdx*m_srcStride, inn...);
  }
  template<class... InN>
  void OddSubsampledAccumulate(I filtIdx, J inIdx, const InN... inn) {
    m_acc+=OddSubsampledAccumulator<T,I,J,Filtn...>::acc(
      filtIdx, inIdx*m_srcStride, inn...);
  }
  void write(I dstIdx) {
    m_dstPtr[dstIdx*m_dstStride]=m_acc;
  }
  void reset() {
    m_acc=0;
  }
  void incrementDstPtr(I dstInc) {
    m_dstPtr+=dstInc;
  }
  void incrementSrcDstPtrStrided(I dstInc) {
    m_dstPtr+=dstInc*m_dstStride;
  }
 
 protected:
  T* m_dstPtr;
  U m_acc;
  I m_srcStride;
  I m_dstStride;
};

template<class Filt, class... Filtn>
struct Updater {
  template<typename I, class Acc, class... AccN>
  static void accumulate(I filtIdx, I srcIdx, Acc&& acc, AccN&& ... accn) {
    Updater<Filt>::accumulate(filtIdx, srcIdx, std::forward<Acc>(acc));
    Updater<Filtn...>::accumulate(filtIdx, srcIdx,std::forward<AccN>(accn)...);
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
  template<typename I, class Acc>
  static void accumulate(I filtIdx, I srcIdx, Acc&& acc) {
    if((filtIdx>=-Filt::TapSizeLeft) && (filtIdx<=Filt::TapSizeRight)) {
      acc.accumulate(srcIdx,Filt::Buff[Filt::TapSizeLeft+filtIdx]);
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

template<typename Acc=void, typename... AccN>
struct SrcDstPtrUpdater {
  static void IncrementSrcDstPtr(size_t srcInc, size_t dstInc,
      Acc&& acc, AccN&&... accn) {
    acc.incrementSrcDstPtr(srcInc, dstInc);
    SrcDstPtrUpdater<AccN...>::IncrementSrcDstPtr(srcInc, dstInc,
      std::forward<AccN>(accn)...);
  }
  static void IncrementSrcDstPtrStrided(size_t srcInc, size_t dstInc,
      Acc&& acc, AccN&&... accn) {
    acc.incrementSrcDstPtrStrided(srcInc, dstInc);
    SrcDstPtrUpdater<AccN...>::IncrementSrcDstPtrStrided(srcInc, dstInc,
      std::forward<AccN>(accn)...);
  }
};
template<>
struct SrcDstPtrUpdater<void> {
  static void IncrementSrcDstPtr(size_t srcInc, size_t dstInc) {}
  static void IncrementSrcDstPtrStrided(size_t srcInc, size_t dstInc) {}
};

template<typename In=void, typename... InN>
struct SimpleUpdater {
  static void Increment(size_t inc, In&& in, InN&&... inn) {
    in += inc;
    SimpleUpdater<InN...>::Increment(inc, std::forward<InN>(inn)...);
  }
};
template<>
struct SimpleUpdater<void> {
  static void Increment(size_t inc) {}
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
  static int PerformSubsampledFilteringXRef(int Nx,
      AccN&&... accn) {

    int Nx_is_odd = (Nx & 1);
    int NxOut = (Nx + Nx_is_odd)/2;

    // Loop over output image
    for (int ox=0; ox<NxOut; ox++) {
      Resetter<AccN...>::reset(std::forward<AccN>(accn)...);
      // Loop over filter size, with periodic boundary conditions
      // TODO TN: this loop can actually be written as a compile time loop
      //#pragma unroll Filt::TapSize
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
        Updater<Filtn...>::accumulate(fx, ix, std::forward<AccN>(accn)...);
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
  template<class Acc, class... InN>
  static int PerformUpsampledFilteringXRef(int NxIn, int NxOut, Acc acc,
    InN... inn) {

    // Loop over output image
    for (int lox=0; lox<NxOut; lox++) {
      int max_x = NxIn-1;
      int ixCentral = lox/2;
      acc.reset();

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
          acc.EvenSubsampledAccumulate(fx, idx_x, inn...);
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
	      acc.OddSubsampledAccumulate(fx, idx_x, inn...);
		}

      }
      // Update each buffer with its respective filter
      acc.write(lox);
    }
    return 1;
  }
};

/** \class SeparableSubsampledConvolutionEngine2D
 * \brief Code for the separable subsampled convolution 2D. This class is a
 * variadic template class, because it can handle multiple filtering for each
 * main loop, assuming the filters have the same size
 *
 * \author Thibault Notargiacomo
 */
//TODO TN: replace second line by first one
//template<typename T, class Filt, class... Filtn>
template<typename T, class... Filtn>
class SeparableSubsampledConvolutionEngine2D {
 public:
  /// Defaulted constructor
  SeparableSubsampledConvolutionEngine2D()=default;
  /// Default destructor
  virtual ~SeparableSubsampledConvolutionEngine2D()=default;

  /// The main method : perform Subsampled convolution on all columns
  template<typename... AccN>
  static int PerformSubsampledFilteringYRef(int Nx, int Ny,
      AccN&&... accn) {
    // We decided to use the tuple trick
    // so that each loop index has its own copy of the accumulator and then
    // can be properly parallelized through openMP for instance
    auto tuple = std::make_tuple(accn...);

    // Loop over output image along y direction
    #pragma omp parallel for
    for (int ox=0; ox<Nx; ox++) {
      // First make a local iteration copy of the accumulator
      std::tie(accn...) = tuple;
      // Second, update the address of buffer
      SrcDstPtrUpdater<AccN...>::IncrementSrcDstPtr(ox, ox,
        std::forward<AccN>(accn)...);
      //Now, you can launch the X convolution
      SeparableSubsampledConvolutionEngine<T, Filtn...
        >::PerformSubsampledFilteringXRef(Ny, std::forward<AccN>(accn)...);
    }
    return 1;
  }

  /// The main method : perform Subsampled convolution on all rows
  template<typename... AccN>
  static int PerformSubsampledFilteringXRef( int Nx, int Ny,
      AccN&&... accn) {
    // We decided to use the tuple trick
    // so that each loop index has its own copy of the accumulator and then
    // can be properly parallelized through openMP for instance
    auto tuple = std::make_tuple(accn...);
 
    // Loop over output image along y direction
    #pragma omp parallel for
    for (int oy=0; oy<Ny; oy++) {
      // First make a local iteration copy of the accumulator
      std::tie(accn...) = tuple;
      // Second, update the address of buffer
      SrcDstPtrUpdater<AccN...>::IncrementSrcDstPtrStrided(oy, oy,
        std::forward<AccN>(accn)...);
      //Now, you can launch the X convolution
      SeparableSubsampledConvolutionEngine<T, Filtn...
        >::PerformSubsampledFilteringXRef(Nx, std::forward<AccN>(accn)...);
    }
    return 1;
  }
};

/** \class SeparableUpsampledConvolutionEngine2D
 * \brief Code for the separable upsampled convolution.
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class... Filtn>
class SeparableUpsampledConvolutionEngine2D {
 public:
  /// Defaulted constructor
  SeparableUpsampledConvolutionEngine2D()=default;
  /// Default destructor
  virtual ~SeparableUpsampledConvolutionEngine2D()=default;

  /// The main method : perform Subsampled convolution on all rows
  template<typename... InN>
  static int PerformUpsampledFilteringXRef( int NxIn, int NxOut,
    int NyIn, int NyOut, T* out, InN&&... inn) {
    // We decided to use the tuple trick
    // so that each loop index has its own copy of the accumulator and then
    // can be properly parallelized through openMP for instance
    auto tuple = std::make_tuple(inn...);

    // Loop over output image along y direction, only X direction will expand
    #pragma omp parallel for
    for (int oy=0; oy<NyIn; oy++) {
      // First make a local iteration copy of the accumulator
      std::tie(inn...) = tuple;
      // Second, update the address of buffer
      SimpleUpdater<InN...>::Increment(oy*NxIn, std::forward<InN>(inn)...);
      //Now, you can launch the X convolution
/*      SeparableUpsampledConvolutionEngine<T,Filtn...
          >::PerformUpsampledFilteringXRef(NxIn, NxOut,
        out+oy*NxOut, inn...);*/
    }
    return 1;
  }
  /// The main method : perform Subsampled convolution on all rows
  template<typename... InN>
  static int PerformUpsampledFilteringYRef( int NxIn, int NxOut,
    int NyIn, int NyOut, T* out, InN&&... inn) {
    // We decided to use the tuple trick
    // so that each loop index has its own copy of the accumulator and then
    // can be properly parallelized through openMP for instance
    auto tuple = std::make_tuple(inn...);

    // Loop over output image along X direction, only Y direction will expand
    #pragma omp parallel for
    for (int ox=0; ox<NxOut; ox++) {
      // First make a local iteration copy of the accumulator
      std::tie(inn...) = tuple;
      // Second, update the address of buffer
      SimpleUpdater<InN...>::Increment(ox, std::forward<InN>(inn)...);
      //Now, you can launch the X convolution
/*      SeparableUpsampledConvolutionEngine<T,Filtn...
          >::PerformUpsampledFilteringXRef(NyIn, NyOut,
        out+ox, inn...);*/
    }
    return 1;
  }};


#endif //SEPARABLE_H
