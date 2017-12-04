// STL
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

// Lib
#include "CoeffContainer.h"
#include "Wavelet3D.h"

//Local
#include "TestEngine.h"

template<
  template<typename> class C,
  typename T,//data type
  typename L,//nb of level
  typename SX,//Size of input data X
  typename SY,//Size of input data Y
  typename SZ>//Size of input data Z
struct Wavelet3DTestFunctor {
  static bool test() {

    //more practical
    constexpr const int level = L::value;
    constexpr const int sizeX = SX::value;
    constexpr const int sizeY = SY::value;
    constexpr const int sizeZ = SZ::value;

    // Define input/output
    std::vector<T> in(sizeX*sizeY*sizeZ);
    // Seed with a real random value, if available
    std::random_device re;
    // Random numbers between -1. and 1.
    std::default_random_engine rnd(re());
    std::uniform_real_distribution<T> uni(-1., 1.);
    auto rand = [&]() { return uni(rnd); };
    std::generate(in.begin(), in.end(),rand);
    const std::vector<T> incopy(in.cbegin(),in.cend());

    // Define wavelet tranform
    C<T> w(in.data(),sizeX,sizeY,sizeZ,false,"Test",level);
    // perform forward transform
    w.forward();
    //Delete previous image
    std::fill(in.begin(),in.end(),0);
    // perform inverse transform
    w.backward(); 
   
    //check reconstruction quality
    return std::inner_product(incopy.cbegin(),incopy.cend(),in.cbegin(), true,
      [](bool acc0, bool acc1) {
        return acc0&&acc1;
      },
      [](T a, T b) {
        return std::abs(a-b) < 1e-4;
      });
    return true;
  }
};

template<typename T, typename L, typename SX, typename SY, typename SZ>
using Daub2_3DTestFctr = Wavelet3DTestFunctor<Daub2_3D,T,L,SX,SY,SZ>;
template<typename T, typename L, typename SX, typename SY, typename SZ>
using Daub3_3DTestFctr = Wavelet3DTestFunctor<Daub3_3D,T,L,SX,SY,SZ>;
template<typename T, typename L, typename SX, typename SY, typename SZ>
using Daub4_3DTestFctr = Wavelet3DTestFunctor<Daub4_3D,T,L,SX,SY,SZ>;
template<typename T, typename L, typename SX, typename SY, typename SZ>
using Daub5_3DTestFctr = Wavelet3DTestFunctor<Daub5_3D,T,L,SX,SY,SZ>;
template<typename T, typename L, typename SX, typename SY, typename SZ>
using Anto97_BiOrth_3DTestFctr = Wavelet3DTestFunctor<
  Anto97_BiOrth_3D,T,L,SX,SY,SZ>;
template<typename T, typename L, typename SX, typename SY, typename SZ>
using QSHIFT6_Orth_3DTestFctr = Wavelet3DTestFunctor<
  QSHIFT6_Orth_3D,T,L,SX,SY,SZ>;
template<typename T, typename L, typename SX, typename SY, typename SZ>
using REVERSE_QSHIFT6_Orth_3DTestFctr = 
  Wavelet3DTestFunctor<REVERSE_QSHIFT6_Orth_3D,T,L,SX,SY,SZ>;
template<typename T, typename L, typename SX, typename SY, typename SZ>
using dtwAnto97QSHIFT6_3DTestFctr =
  Wavelet3DTestFunctor<dtwAnto97QSHIFT6_3D,T,L,SX,SY,SZ>;  


int main(int argc, char* argv[])  {

  // Defining different sizes
  using T = std::tuple<float,double>;
  using L = std::tuple<IntType<1>,IntType<2>,IntType<3>,IntType<4>>;
  using S = std::tuple<IntType<63>,IntType<63>,IntType<64>,IntType<65>,
    IntType<66>, IntType<127>>;

  //We challenge the template test functor over the
  //cartesian product of the type sets T,L,S
  return 
    TestImp<Daub2_3DTestFctr>(T(), L(), S(), S(), S()) &&
    TestImp<Daub3_3DTestFctr>(T(), L(), S(), S(), S()) &&
    TestImp<Daub4_3DTestFctr>(T(), L(), S(), S(), S()) &&
    TestImp<Daub5_3DTestFctr>(T(), L(), S(), S(), S()) &&
    TestImp<Anto97_BiOrth_3DTestFctr>(T(), L(), S(), S(), S()) &&
    TestImp<QSHIFT6_Orth_3DTestFctr>(T(), L(), S(), S(), S()) &&
    TestImp<REVERSE_QSHIFT6_Orth_3DTestFctr>(T(), L(), S(), S(), S()) &&
    TestImp<dtwAnto97QSHIFT6_3DTestFctr>(T(), L(), S(), S(), S()) 
    ? EXIT_SUCCESS : EXIT_FAILURE;
    EXIT_SUCCESS;
}
