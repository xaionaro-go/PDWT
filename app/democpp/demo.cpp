// STL
#include <cstdlib>
#include <numeric>
#include <vector>

// Lib
#include "coeffContainer.h"
#include "Wavelet1D.h"

int main(int argc, char **argv) {

  // Define input
  std::vector<float> in(10);
  std::iota(in.begin(), in.end(),0);

  CoeffContainer1D<float,std::vector<float,PackAllocator<float>>> a;
  // Define wavelet tranform
  //Daub2_1D<float> w(
  //  in.data(),in.size(),1,1,false,"Daub2",3);
  //w.forward();
  //w.backward();
  return EXIT_SUCCESS;
}
