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
      //std::cout<<" , "<<2*(filtIdx+Filt::TapHalfSizeLeft)+
      //  Filt::EvenSubSampOffset<<")="<<inIdx;
	  return Filt::Buff[2*(filtIdx+Filt::TapHalfSizeLeft)+
        Filt::EvenSubSampOffset]*in[inIdx];
    } else {
      return (T)0;
    }
  }
};
template<typename T, typename I, typename J>
struct EvenSubsampledAccumulator<T,I,J,void> {
 static constexpr T acc(I filtIdx, J inIdx) { return (T)0; }
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
      //std::cout<<2*(filtIdx+Filt::TapHalfFloorSizeLeft)+
      //  Filt::OddSubSampOffset<<")="<<inIdx;
	  return Filt::Buff[2*(filtIdx+Filt::TapHalfFloorSizeLeft)+
        Filt::OddSubSampOffset]*in[inIdx];
    } else {
      return (T)0;
    }
  }
};
template<typename T, typename I, typename J>
struct OddSubsampledAccumulator<T,I,J,void> {
  static constexpr T acc(I filtIdx, J inIdx) { return (T)0; }
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
  void EvenSubsampledAccumulate(I filtIdx, J inIdx, InN&&... inn) {
    m_acc+=EvenSubsampledAccumulator<T,I,J,Filtn...>::acc(
      filtIdx, inIdx*m_srcStride, std::forward<InN>(inn)...);
  }
  template<class... InN>
  void OddSubsampledAccumulate(I filtIdx, J inIdx, InN&&... inn) {
    m_acc+=OddSubsampledAccumulator<T,I,J,Filtn...>::acc(
      filtIdx, inIdx*m_srcStride, std::forward<InN>(inn)...);
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

template<typename T, typename U, typename I, typename J, class... Filtn>
class SubsampledAccumulatorUpdate : 
    public SubsampledAccumulator<T,U,I,J,Filtn...> {
 public:
  SubsampledAccumulatorUpdate(T* dstPtr, I srcStride=1, I dstStride=1):
    SubsampledAccumulator<T,U,I,J,Filtn...>(dstPtr , srcStride, dstStride)  {}

  void write(I dstIdx) {
    this->m_dstPtr[dstIdx*this->m_dstStride]+=this->m_acc;
  }
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
//TODO TN      std::cout<<"Filter is "<<Filt::Buff[Filt::TapSizeLeft+filtIdx]<<std::endl;
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

/** \class SeparableSubsampledConvolutionEngine
 * \brief Code for the separable subsampled convolution. This class is a
 * variadic template class, because it can handle multiple filtering for each
 * main loop, assuming the filters have the same size
 *
 * TODO TN: perf issue: you should use temporary for accumulation
 *
 * \author Thibault Notargiacomo
 */
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

      //std::cout<<"Lox="<<lox<<" : ";
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
          //std::cout<<"("<<idx_x;
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
          //std::cout<<"("<<idx_x<<" , ";
	      acc.OddSubsampledAccumulate(fx, idx_x, inn...);
		}
      }
      //std::cout<<std::endl;
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
template<typename T, class... Filtn>
class SeparableSubsampledConvolutionEngine2D {
 public:
  /// Defaulted constructor
  SeparableSubsampledConvolutionEngine2D()=default;
  /// Default destructor
  virtual ~SeparableSubsampledConvolutionEngine2D()=default;

  /// The main method : perform Subsampled convolution on all columns
  template<typename... AccN>
  static int PerformSubsampledFilteringYRef(int NxIn, int NyIn,
      AccN&&... accn) {
    // Loop over output image along y direction
    #pragma omp parallel for
    for (int ox=0; ox<NxIn; ox++) {
      callSubsampledConvWithTuple(
        NyIn,
        ox,
        ox,
        std::make_tuple(accn...),
        std::index_sequence_for<AccN...>());
    }
    return 1;
  }

  /// The main method : perform Subsampled convolution on all rows
  template<typename... AccN>
  static int PerformSubsampledFilteringXRef( int NxIn, int NxOut, int Ny,
      AccN&&... accn) {
    // Loop over output image along y direction
    #pragma omp parallel for
    for (int oy=0; oy<Ny; oy++) {
      //Now, you can launch the X convolution
      callSubsampledConvWithTuple(
        NxIn,
        oy*NxIn,
        oy*NxOut,
        std::make_tuple(accn...),
        std::index_sequence_for<AccN...>());
    }
    return 1;
  }
 protected:
  template<typename... AccN, std::size_t... Is>
  static void callSubsampledConvWithTuple(
      int Nin, int srcStride, int dstStride,
      std::tuple<AccN...>&& accn, std::index_sequence<Is...>) {
    SrcDstPtrUpdater<AccN...>::IncrementSrcDstPtr(
      srcStride, dstStride, std::get<Is>(std::forward<std::tuple<AccN...>>(accn))...);
    SeparableSubsampledConvolutionEngine<T, Filtn...
        >::PerformSubsampledFilteringXRef(Nin, std::get<Is>(accn)...);
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
      int Ny, T* out, InN... inn) {
    // Loop over output image along y direction, only X direction will expand
    #pragma omp parallel for
    for (int oy=0; oy<Ny; oy++) {
      //Now, you can launch the X convolution
      callUpsampledConvWithTuple(
        NxIn,
        NxOut,
        SubsampledAccumulator<T,T,int,int,Filtn...>(out+oy*NxOut),
        std::make_tuple((inn+oy*NxIn)...),
        std::index_sequence_for<InN...>());
    }
    return 1;
  }
  /// The main method : perform Subsampled convolution on all rows
  template<typename... InN>
  static int PerformUpsampledFilteringYRef( int Nx, int NyIn, int NyOut,
      T* out, InN... inn) {
    #pragma omp parallel for
    for (int ox=0; ox<Nx; ox++) {
      callUpsampledConvWithTuple(
        NyIn,
        NyOut,
        SubsampledAccumulator<T,T,int,int,Filtn...>(out+ox,Nx,Nx),
        std::make_tuple((inn+ox)...),
        std::index_sequence_for<InN...>());
    }
    return 1;
  }
 protected:
  template<typename Acc, typename... Inputs, std::size_t... Is>
  static void callUpsampledConvWithTuple(int Nin, int Nout, Acc&& acc,
      std::tuple<Inputs...>&& inn, std::index_sequence<Is...>) {
    SeparableUpsampledConvolutionEngine<T,Filtn...
        >::PerformUpsampledFilteringXRef(
      Nin,
      Nout,
      acc,
      std::get<Is>(inn)...);
  }
};

/** \class SeparableSubsampledConvolutionEngine3D
 * \brief Code for the separable subsampled convolution 3D. This class is a
 * variadic template class, because it can handle multiple filtering for each
 * main loop, assuming the filters have the same size
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class... Filtn>
class SeparableSubsampledConvolutionEngine3D {
 public:
  /// Defaulted constructor
  SeparableSubsampledConvolutionEngine3D()=default;
  /// Default destructor
  virtual ~SeparableSubsampledConvolutionEngine3D()=default;

  /// The main method : perform Subsampled convolution on all columns
  template<typename... AccN>
  static int PerformSubsampledFilteringZRef(int NxIn, int NyIn, int NzIn,
      AccN&&... accn) {
    size_t yStride=NxIn;
    // Loop over both y and x to perform z filtering
    for (int oy=0; oy<NyIn; oy++) {
      #pragma omp parallel for
      for (int ox=0; ox<NxIn; ox++) {
        callSubsampledConvWithTuple(
          NzIn,
          oy*yStride+ox,
          oy*yStride+ox,
          std::make_tuple(accn...),
          std::index_sequence_for<AccN...>());
      }
    }
    return 1;
  }

  /// The main method : perform Subsampled convolution on all columns
  template<typename... AccN>
  static int PerformSubsampledFilteringYRef(int NxIn, int NyIn, int NyOut,
    int NzIn, AccN&&... accn) {
    size_t zStrideIn=NyIn*NxIn;
    size_t zStrideOut=NyOut*NxIn;
    // Loop over both z and x to perform y filtering
    for (int oz=0; oz<NzIn; oz++) {
      #pragma omp parallel for
      for (int ox=0; ox<NxIn; ox++) {
        callSubsampledConvWithTuple(
          NyIn,
          oz*zStrideIn+ox,
          oz*zStrideOut+ox,
          std::make_tuple(accn...),
          std::index_sequence_for<AccN...>());
      }
    }
    return 1;
  }

  /// The main method : perform Subsampled convolution on all rows
  template<typename... AccN>
  static int PerformSubsampledFilteringXRef( int NxIn, int NxOut,
    int NyIn, int NzIn, AccN&&... accn) {
    size_t zStrideIn=NxIn*NyIn;
    size_t zStrideOut=NxOut*NyIn;
    // Loop over both z and y to perform x filtering
    for (int oz=0; oz<NzIn; oz++) {
      #pragma omp parallel for
      for (int oy=0; oy<NyIn; oy++) {
        //Now, you can launch the X convolution
        callSubsampledConvWithTuple(
          NxIn,
          oz*zStrideIn+oy*NxIn,
          oz*zStrideOut+oy*NxOut,
          std::make_tuple(accn...),
          std::index_sequence_for<AccN...>());
      }
    }
    return 1;
  }
 protected:
  template<typename... AccN, std::size_t... Is>
  static void callSubsampledConvWithTuple(
      int Nin, int srcStride, int dstStride,
      std::tuple<AccN...>&& accn, std::index_sequence<Is...>) {
    SrcDstPtrUpdater<AccN...>::IncrementSrcDstPtr(
      srcStride, dstStride, std::get<Is>(std::forward<std::tuple<AccN...>>(accn))...);
    SeparableSubsampledConvolutionEngine<T, Filtn...
        >::PerformSubsampledFilteringXRef(Nin, std::get<Is>(accn)...);
  }
};

/** \class SeparableUpsampledConvolutionEngine3D
 * \brief Code for the separable upsampled convolution.
 *
 * \author Thibault Notargiacomo
 */
template<typename T,
  template <typename,typename,typename,typename, class...> class AccT,
  class... Filtn>
class SeparableUpsampledConvolutionEngine3D {
 public:
  /// Defaulted constructor
  SeparableUpsampledConvolutionEngine3D()=default;
  /// Default destructor
  virtual ~SeparableUpsampledConvolutionEngine3D()=default;

  /// The main method : perform Upsampled convolution on all rows
  template<typename... InN>
  static int PerformUpsampledFilteringXRef( int NxIn, int NxOut,
      int Ny, int Nz, T* out, InN... inn) {
    size_t zStrideIn=Ny*NxIn;
    size_t zStrideOut=Ny*NxOut;
    // Loop over output image along y direction, only X direction will expand
    for (int oz=0; oz<Nz; oz++) {
      #pragma omp parallel for
      for (int oy=0; oy<Ny; oy++) {
        //Now, you can launch the X convolution
        callUpsampledConvWithTuple(
          NxIn,
          NxOut,
          AccT<T,T,int,int,Filtn...>(out+oz*zStrideOut+oy*NxOut),
          std::make_tuple((inn+oz*zStrideIn+oy*NxIn)...),
          std::index_sequence_for<InN...>());
      }
    }
    return 1;
  }
  /// The main method : perform Upsampled convolution on all rows
  template<typename... InN>
  static int PerformUpsampledFilteringYRef( int Nx, int NyIn, int NyOut,
      int Nz, T* out, InN... inn) {
    size_t zStrideIn=NyIn*Nx;
    size_t zStrideOut=NyOut*Nx;
    for (int oz=0; oz<Nz; oz++) {
      #pragma omp parallel for
      for (int ox=0; ox<Nx; ox++) {
        callUpsampledConvWithTuple(
          NyIn,
          NyOut,
          AccT<T,T,int,int,Filtn...>(out+oz*zStrideOut+ox,Nx,Nx),
          std::make_tuple((inn+oz*zStrideIn+ox)...),
          std::index_sequence_for<InN...>());
      }
    }
    return 1;
  }
  /// Upsampled convolution for z direction
  template<typename... InN>
  static int PerformUpsampledFilteringZRef( int Nx, int Ny,
      int NzIn, int NzOut, T* out, InN... inn) {
    size_t zStride=Ny*Nx;
    for (int oy=0; oy<Ny; oy++) {
      #pragma omp parallel for
      for (int ox=0; ox<Nx; ox++) {
        callUpsampledConvWithTuple(
          NzIn,
          NzOut,
          AccT<T,T,int,int,Filtn...>(out+oy*Nx+ox,zStride,zStride),
          std::make_tuple((inn+oy*Nx+ox)...),
          std::index_sequence_for<InN...>());
      }
    }
    return 1;
  }

 protected:
  template<typename Acc, typename... Inputs, std::size_t... Is>
  static void callUpsampledConvWithTuple(int Nin, int Nout, Acc&& acc,
      std::tuple<Inputs...>&& inn, std::index_sequence<Is...>) {
    SeparableUpsampledConvolutionEngine<T,Filtn...
        >::PerformUpsampledFilteringXRef(
      Nin,
      Nout,
      acc,
      std::get<Is>(inn)...);
  }
};

#endif //SEPARABLE_H
