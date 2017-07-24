// STL
#include <cmath>


/** \struct WaveletToCpx2D
 * \brief Equivalent of a matrix multiplication for complex interpretation
 *
 * This functor helps to convert from the 4 output vectors
 * of dtcwt filtering to 4 vectors that can be considered as
 * 2 complex vectors of 2 coefficients (real - imaginary).
 *
 * This operator is in fact equivalent to the following scaled orthogonal
 * matrix : <br>
 *
 *	\f$ \frac{1}{\sqrt{2}} \times
		\begin{pmatrix}  1 & 0 & 0 & 1 \\
 	 	 	 	 	 	 0 & 1 & 1 & 0 \\
 	 	 	 	 	 	 0 & 1 & -1 & 0 \\
 	 	 	 	 	 	 1 & 0 & 0 & -1 \end{pmatrix} \f$
 *
 * \author Thibault Notargiacomo
 */
template<typename T>
struct WaveletToCpx2D {
  WaveletToCpx2D( T ratio ) : m_ratio( ratio ) {}
  template <class Tuple>
  void operator()(Tuple in) const {
    T tmp = (boost::get<0>(in)+boost::get<3>(in)) * m_ratio;
    boost::get<3>(in) = (boost::get<0>(in)-boost::get<3>(in))*m_ratio;
    boost::get<0>(in) = tmp;
    tmp               = (boost::get<1>(in)+boost::get<2>(in))*m_ratio;
    boost::get<2>(in) = (boost::get<1>(in)-boost::get<2>(in))*m_ratio;
    boost::get<1>(in) = tmp;
  }
  const T m_ratio;
};

/** \struct WaveletToCpx3D
 * \brief Equivalent of a matrix multiplication for complex interpretation
 *
 * This functor helps to convert from the 8 output vectors
 * of dtcwt filtering to 8 vectors that can be considered as
 * 4 complex vectors of 2 coefficients (real - imaginary).
 *
 * This operator is in fact equivalent to the following scaled orthogonal matrix : <br>
 *
 *	\f$ \frac{1}{4\sqrt{2}} \times \begin{pmatrix}
       1 & 0 &  0 & -1 &  0 & -1 & -1 &  0 \\
       0 & 1 &  1 &  0 &  1 &  0 &  0 & -1 \\
       1 & 0 &  0 & -1 &  0 &  1 &  1 &  0 \\
       0 & 1 &  1 &  0 & -1 &  0 &  0 &  1 \\
       1 & 0 &  0 &  1 &  0 & -1 &  1 &  0 \\
       0 & 1 & -1 &  0 &  1 &  0 &  0 &  1 \\
       1 & 0 &  0 &  1 &  0 &  1 & -1 &  0 \\
       0 & 1 & -1 &  0 & -1 &  0 &  0 & -1 \end{pmatrix} \f$
 *
 * \author Thibault Notargiacomo
 *
template<typename T>
struct WaveletToCpx3D {
  WaveletToCpx3D( T ratio ) : m_ratio( ratio ) {}

  template <class Tuple>
  void operator()(Tuple in) const {
    T T0 = boost::get<0>(in);
    T T1 = boost::get<1>(in);
    T T2 = boost::get<2>(in);
    T T3 = boost::get<3>(in);
    T T4 = boost::get<4>(in);
    T T5 = boost::get<5>(in);
    T T6 = boost::get<6>(in);
    T T7 = boost::get<7>(in);

    //Treating Real Part
    boost::get<0>(in) 	= (T0-T3-T5-T6) * m_ratio;
    boost::get<2>(in) 	= (T0-T3+T5+T6) * m_ratio;
    boost::get<4>(in) 	= (T0+T3-T5+T6) * m_ratio;
    boost::get<6>(in) 	= (T0+T3+T5-T6) * m_ratio;

    //Treating Imaginary Part
    boost::get<1>(in) 	= (T1+T2+T4-T7) * m_ratio;
    boost::get<3>(in) 	= (T1+T2-T4+T7) * m_ratio;
    boost::get<5>(in) 	= (T1-T2+T4+T7) * m_ratio;
    boost::get<7>(in) 	= (T1-T2-T4-T7) * m_ratio;
  }
  const T m_ratio;
};


/** \struct CpxToWavelet3D
 * \brief Equivalent of a matrix multiplication for complex interpretation
 *
 * This functor helps to convert from the 8
 * complex vectors of coefficients in the dtcwt basis
 * to 8 vectors of coefficients that will be ready
 * for the synthesis operation.
 *
 * This operator is in fact equivalent to the following
 * scaled orthogonal matrix : <br>
 *
 *	\f$ \frac{1}{4\sqrt{2}} \times
    \begin{pmatrix}
      1 &  0 &  1 &  0 &  1 &  0 &  1 &  0 \\
      0 &  1 &  0 &  1 &  0 &  1 &  0 &  1 \\
      0 &  1 &  0 &  1 &  0 & -1 &  0 & -1 \\
      -1 &  0 & -1 &  0 &  1 &  0 &  1 &  0 \\
      0 &  1 &  0 & -1 &  0 &  1 &  0 & -1 \\
      -1 &  0 &  1 &  0 & -1 &  0 &  1 &  0 \\
      -1 &  0 &  1 &  0 &  1 &  0 & -1 &  0 \\
      0 & -1 &  0 &  1 &  0 &  1 &  0 & -1 \end{pmatrix} \f$
 *
 * \author Thibault Notargiacomo
 *
template<typename T>
struct CpxToWavelet3D {
  CpxToWavelet3D( T ratio ) : m_ratio( ratio ) {}

  template <class Tuple>
  void operator()(Tuple in) const {
    T T0 = boost::get<0>(in);
    T T1 = boost::get<1>(in);
    T T2 = boost::get<2>(in);
    T T3 = boost::get<3>(in);
    T T4 = boost::get<4>(in);
    T T5 = boost::get<5>(in);
    T T6 = boost::get<6>(in);
    T T7 = boost::get<7>(in);

    //Treating Real Part
    boost::get<0>(in) 	= (T0+T2+T4+T6) * m_ratio;
    boost::get<1>(in) 	= (T1+T3+T5+T7) * m_ratio;
    boost::get<2>(in) 	= (T1+T3-T5-T7) * m_ratio;
    boost::get<3>(in) 	= (-T0-T2+T4+T6) * m_ratio;
    boost::get<4>(in) 	= (T1-T3+T5-T7) * m_ratio;
    boost::get<5>(in) 	= (-T0+T2-T4+T6) * m_ratio;
    boost::get<6>(in) 	= (-T0+T2+T4-T6) * m_ratio;
    boost::get<7>(in) 	= (-T1+T3+T5-T7) * m_ratio;
  }
  const T m_ratio;
};
*/
