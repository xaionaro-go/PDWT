// STL
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

// Lib
#include "coeffContainer.h"
#include "Wavelet1D.h"

int main(int argc, char **argv) {

  auto print = [](auto& in){ std::cout<<in<<", "; };

  // Define input
  std::vector<float> in(10);
  std::iota(in.begin(), in.end(),0);

  std::cout<<"Input is: ";
  std::for_each(in.begin(),in.end(),print);
  std::cout<<std::endl;

  // Define wavelet tranform
  Daub2_1D<float> w(in.data(),in.size(),1,1,false,"Daub2",3);
  w.forward();
 
  std::cout<<"Output is: ";
  std::for_each(w.get_coeff().begin(),w.get_coeff().end(),print);
  std::cout<<std::endl;

  w.backward();
  return EXIT_SUCCESS;
}
