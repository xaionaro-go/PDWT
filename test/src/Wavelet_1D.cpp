// STL
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <vector>

// Lib
#include "coeffContainer.h"
#include "Wavelet1D.h"

//Local
#include "TestEngine.h"

template<
  template<typename> class C,
  typename T,//data type
  typename L,//nb of level
  typename S>//Size of input data
struct Wavelet1DTestFunctor {
  static bool test() {

    //more practical
    constexpr const int level = L::value;
    constexpr const int size = S::value;

    // Define input/output
    std::vector<T> in(size);
    std::iota(in.begin(), in.end(),0);
    const std::vector<T> incopy(in.cbegin(),in.cend());
    std::vector<T> out(in.size());

    // Define wavelet tranform
    C<T> w(in.data(),in.size(),1,1,false,"Test",level);
    // perform forward transform
    w.forward();
    // perform inverse transform
    w.backward(); 
    //get reconstruction 
    w.get_image(out.data());
   
    //check reconstruction quality
    return std::inner_product(incopy.cbegin(),incopy.cend(),out.cbegin(), true,
      [](bool acc0, bool acc1) {
        return acc0&&acc1;
      },
      [](T a, T b) {
        return std::abs(a-b) < 1;
      });
  }
};

template<typename T, typename L, typename S>
using Daub2_1DTestFctr = Wavelet1DTestFunctor<Daub2_1D,T,L,S>;
template<typename T, typename L, typename S>
using Daub3_1DTestFctr = Wavelet1DTestFunctor<Daub3_1D,T,L,S>;
template<typename T, typename L, typename S>
using Daub4_1DTestFctr = Wavelet1DTestFunctor<Daub4_1D,T,L,S>;
template<typename T, typename L, typename S>
using Daub5_1DTestFctr = Wavelet1DTestFunctor<Daub5_1D,T,L,S>;
template<typename T, typename L, typename S>
using Anto97_BiOrth_1DTestFctr = Wavelet1DTestFunctor<Anto97_BiOrth_1D,T,L,S>;
template<typename T, typename L, typename S>
using QSHIFT6_Orth_1DTestFctr = Wavelet1DTestFunctor<QSHIFT6_Orth_1D,T,L,S>;
template<typename T, typename L, typename S>
using REVERSE_QSHIFT6_Orth_1DTestFctr = 
  Wavelet1DTestFunctor<REVERSE_QSHIFT6_Orth_1D,T,L,S>;
template<typename T, typename L, typename S>
using dtwAnto97QSHIFT6_1DTestFctr =
  Wavelet1DTestFunctor<dtwAnto97QSHIFT6_1D,T,L,S>;  

int main(int argc, char* argv[])  {

  // Defining different sizes
  using T = std::tuple<float,double>;
  using L = std::tuple<IntType<1>,IntType<2>,IntType<3>,IntType<4>>;
  using S = std::tuple<IntType<63>,IntType<63>,IntType<64>,IntType<65>,
    IntType<66>, IntType<127>>;

  //We challenge the template test functor over the
  //cartesian product of the type sets T,L,S
  return 
    TestImp<Daub2_1DTestFctr>(T(), L(), S()) &&
    TestImp<Daub3_1DTestFctr>(T(), L(), S()) &&
    TestImp<Daub4_1DTestFctr>(T(), L(), S()) &&
    TestImp<Daub5_1DTestFctr>(T(), L(), S()) &&
    TestImp<Anto97_BiOrth_1DTestFctr>(T(), L(), S()) &&
    TestImp<QSHIFT6_Orth_1DTestFctr>(T(), L(), S()) &&
    TestImp<REVERSE_QSHIFT6_Orth_1DTestFctr>(T(), L(), S()) &&
    TestImp<dtwAnto97QSHIFT6_1DTestFctr>(T(), L(), S()) 
    ? EXIT_SUCCESS : EXIT_FAILURE;
}
