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

  // Define wavelet tranform
  Daub2_1D<float> w(in.data(),in.size(),1,1,false,"Daub2",3);
  w.forward();
  w.backward();

  std::vector<float> out(in.size(), 0);
  w.get_image(out.data());

  //std::all_of(c.cbegin(),c.cend(),[val](float in) { return in==val; });
  bool isOK = std::equal(in.cbegin(),in.cend(),out.cbegin());
  return isOK ? EXIT_SUCCESS:EXIT_FAILURE;
}
