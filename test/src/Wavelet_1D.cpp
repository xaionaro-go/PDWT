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

/*template<
  typename T,
  typename M,
  typename N,
  typename K,
  typename D>
struct Wavelet1DTestFunctor {
  static bool test() {

    //more practical
    constexpr const int m = M::value;
    constexpr const int n = N::value;
    constexpr const int k = K::value;

    // Define input/output
    std::vector<T> in();
    std::vector<T> out(in.size());
    std::iota(in.begin(), in.end(),0);

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
    //std::all_of(c.cbegin(),c.cend(),[val](float in) { return in==val; });

    //Allocate the A, B, C and Ccheck matrices
    std::vector<T> A(m*k);
    std::vector<T> B(k*n);
    std::vector<T> C(m*n);
    std::vector<T> Ccheck(C.size(),0);

    //Init A and B
    Initializer<T> init;
    init.initialize(m, k, n, A.data(), B.data(), C.data());

    //Compute naive serial MM
    for(size_t i=0;i<m;i++){
      for(size_t j=0;j<n;j++){
        for(size_t l=0;l<k;l++){
          Ccheck[j*m+i]+=A[l*m+i]*B[j*k+l];
        }
      }
    }

    //Custom method to be tested
    Carma::Multiply<T>(m, k, n,
      A.data(), B.data(), C.data(), D::value);

    //Correctness checking
    double tol = 1e4*std::numeric_limits<T>::epsilon();
    double refNrm2 = std::sqrt(std::inner_product(Ccheck.cbegin(),
      Ccheck.cend(), Ccheck.cbegin(), 0.0));
    assert(("Norm of the resulting matrix C should not be too small for the\
      correctness test", refNrm2 > std::sqrt(tol)));
    //Now substract Ccheck from C
    std::transform(C.begin(),C.end(),Ccheck.cbegin(),C.begin(),
      std::minus<T>());
    double epsNrm2 = std::sqrt(std::inner_product(C.cbegin(), C.cend(),
      C.cbegin(), 0.0));

    if( epsNrm2/refNrm2 > tol ) {
      std::cout<< "FAILURE: error in matrix multiply exceeds an acceptable "
        <<"margin (m,n,k)=("<<m<<","<<n<<","<<k<<")"<<std::endl;
      return false;
    }
 
    return true;
  }
};
*/

int main(int argc, char* argv[])  {

  // Defining different sizes
  using M = std::tuple<IntType<1>,IntType<2>,IntType<3>,IntType<4>,IntType<5>,
    IntType<64>>;
  using N = std::tuple<IntType<1>,IntType<2>,IntType<3>,IntType<4>,IntType<5>,
    IntType<64>>;
  using K = std::tuple<IntType<1>,IntType<2>,IntType<3>,IntType<4>,IntType<5>,
    IntType<64>>;

  // Defining max recursion size
  using D = std::tuple<IntType<0>,IntType<1>,IntType<2>,
    IntType<3>,IntType<100>>;

  // Carma only defined for single and doubles for now
  using T = std::tuple<float,double>;

  //We challenge the template test functor over the
  //cartesian product of the type sets M,N,K,D
  //return TestImp<MMTestFunctor>(T(), M(), N(), K(), D()) ?
  //  EXIT_SUCCESS : EXIT_FAILURE;
  return EXIT_SUCCESS;
}
