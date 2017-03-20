
// Local
#include "coeffContainer.h"

// STL
#include <vector>

// Local
#include "cuda/managedAllocator.cu.h"


// Explicit instanciation
template class
CoeffContainer1D<float,std::vector<float,managedAllocator<float>>>;
template class
CoeffContainer2D<float,std::vector<float,managedAllocator<float>>>;
template class
CoeffContainer3D<float,std::vector<float,managedAllocator<float>>>;

template class
CoeffContainer1D<double,std::vector<double,managedAllocator<double>>>;
template class
CoeffContainer2D<double,std::vector<double,managedAllocator<double>>>;
template class
CoeffContainer3D<double,std::vector<double,managedAllocator<double>>>;

