// STL
#include <cstdlib>
#include <limits>
#include <numeric>
#include <vector>

// Lib
#include "coeffContainer.h"
#include "Wavelet1D.h"

//Local
#include "TestEngine.h"

/*  auto print = [](auto& in){ std::cout<<in<<", "; };

  // Define input/output
  std::vector<float> in(17);
  std::vector<float> out(in.size());
  std::iota(in.begin(), in.end(),0);
  //std::fill(in.begin(), in.end(), 2);

  std::cout<<"Input is: ";
  std::for_each(in.begin(),in.end(),print);
  std::cout<<std::endl;

  // Define wavelet tranform
  Daub3_1D<float> w(in.data(),in.size(),1,1,false,"Daub2",2);
  // print coeffs when initialized
  std::for_each(w.get_coeff().begin(),w.get_coeff().end(),print);
  // perform forward transform
  w.forward();
  //Delete previous image
  std::fill(in.begin(),in.end(),0);
  // print coeffs
  std::for_each(w.get_coeff().begin(),w.get_coeff().end(),print);
  // perform inverse transform
  w.backward(); 
  
  w.get_image(out.data());
  std::cout<<"Output is: ";
  std::for_each(out.cbegin(),out.cend(),print);
*/

int main(int argc, char **argv) {

  //std::all_of(c.cbegin(),c.cend(),[val](float in) { return in==val; });
  return true?EXIT_SUCCESS:EXIT_FAILURE;
}
