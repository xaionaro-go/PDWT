
// Local
#include "CoeffContainer.h"

// STL
#include <vector>

// Local
#include "vectorization/Vectorization.h"


// Explicit instanciation
template class
CoeffContainer1D<float,std::vector<float,PackAllocator<float>>>;
template class
CoeffContainer2D<float,std::vector<float,PackAllocator<float>>>;
template class
CoeffContainer3D<float,std::vector<float,PackAllocator<float>>>;

template class
CoeffContainer1D<double,std::vector<double,PackAllocator<double>>>;
template class
CoeffContainer2D<double,std::vector<double,PackAllocator<double>>>;
template class
CoeffContainer3D<double,std::vector<double,PackAllocator<double>>>;

