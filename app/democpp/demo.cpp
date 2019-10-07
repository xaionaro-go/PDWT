// STL
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

// Lib
#include "CoeffContainer.h"
#include "Wavelet1D.h"

using T = float;

int main(int argc, char **argv) {

  auto print = [](auto& in){ std::cout<<in<<", "; };

  // Define input/output
  std::vector<T> in(8);
  std::iota(in.begin(), in.end(),0);
  //std::fill(in.begin(), in.end(), 0);

  std::cout<<"Input is: ";
  std::for_each(in.begin(),in.end(),print);
  std::cout<<std::endl;

  // Define wavelet tranform
  Daub2_1D<T> w(in.data(),in.size(),1,1,false,"Daub2",1);
  //Anto97_BiOrth_1D<T> w(in.data(),in.size(),1,1,false,"Anto97",1);
  //REVERSE_QSHIFT6_Orth_1D<T> w(in.data(),in.size(),1,1,false,"QSHIFT6",1);
  //dtwAnto97QSHIFT6_1D<T> w(in.data(),in.size(),1,1,false,"DTCWT",3);  
  //Dummy2_1D<T> w(in.data(),in.size(),1,1,false,"Daub2",1);
 
  // print coeffs when initialized
  std::cout<<"Coefficient after initialization (should be 0)"<<std::endl;
  std::for_each(w.get_coeff().begin(),w.get_coeff().end(),print);
  std::cout<<std::endl;

  // perform forward transform
  w.forward();
  //Delete previous image
  std::fill(in.begin(),in.end(),0);

  // print coeffs
  std::cout<<"Coefficient after transform (should NOT be 0)"<<std::endl;
  std::for_each(w.get_coeff().begin(),w.get_coeff().end(),print);
  std::cout<<std::endl;

  // perform inverse transform
  w.backward(); 
  
  std::cout<<"Output should now be equal to previous input: ";
  std::for_each(in.cbegin(),in.cend(),print);


  std::cout<<std::endl;

  return EXIT_SUCCESS;
}
