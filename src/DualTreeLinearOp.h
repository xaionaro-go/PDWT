// STL
#include <cmath>


/** \struct WaveletToCpx2D
 * \brief A linear operator that allows to perform the alternative x/y negative
 * frequency cancellation
 *
 * From the 4 raw output of the dtcwt filtering process, one actually obtains
 * a combination of signals(real): a , and their hilbert transform(imag): b
 * We recall that the complex resulting vector 1/2 (a + i b) has its negative
 * frequency content cancelled.
 * 1/2 (a - ib) has its positive frequency cancelled
 *
 * More precisely, let say we have 4 row vectors as input, they
 * would represent:
 * 0: real Y * real X : real
 * 1: real Y * Imag X : imag
 * 2: imag Y * real X : imag
 * 3: imag Y * Imag X : -real
 *
 * When doing a combination of real and Hilbert/imaginary filters along x and y
 * where negative frequency cancellation property is used, we obtain a
 * combination that has a simple but limited interpretation:
 * (realY+iImagY)(realX+iImagX)
 * =realY realX - ImagY ImagX + i (realY ImagX + ImagY RealX)
 * =(0-3) + i(1+2) : positive X and positive Y freq quadrant
 *
 * One can use a negative y cancellation, followed by a x positive cancellation
 * (realY-iImagY)*(realX+iImagX):
 * = realY realX + ImagY ImagX + i (RealY ImagX - ImagY RealX)
 * (0+3) + i(1-2)
 *
 * Resulting operator can be written: <br>
 *
 *	\f$ \frac{1}{\sqrt{2}} \times
        \begin{pmatrix}  1 & 0 & 0 & -1 \\
                         0 & 1 & 1 & 0 \\
                         1 & 0 & 0 & 1 \\
                         0 & 1 & -1 & 0 \\
                         \end{pmatrix} \f$
 *
 * \author Thibault Notargiacomo
 */
template<typename T>
struct WaveletToCpx2D {
  WaveletToCpx2D( T ratio ) : m_ratio( ratio ) {}
  template <class Tuple>
  void operator()(Tuple in) const {
    T tmp2 = boost::get<2>(in);
    boost::get<2>(in) = (boost::get<0>(in)+boost::get<3>(in))*m_ratio;
    boost::get<0>(in) = (boost::get<0>(in)-boost::get<3>(in))*m_ratio;
    boost::get<3>(in) = (boost::get<1>(in)-tmp2)*m_ratio;
    boost::get<1>(in) = (boost::get<1>(in)+tmp2)*m_ratio;
  }
  const T m_ratio;
};

/** \struct CpxToWavelet2D
 * \brief Just the inverse of WaveletToCpx2D
 *
 * Resulting operator can be written: <br>
 *
 *	\f$ \frac{1}{\sqrt{2}} \times
        \begin{pmatrix}  1 & 0 & 1 & 0 \\
                         0 & 1 & 0 & 1 \\
                         0 & 1 & 0 & -1 \\
                         -1 & 0 & 1 & 0 \\
                         \end{pmatrix} \f$
 *
 * \author Thibault Notargiacomo
 */
template<typename T>
struct CpxToWavelet2D {
  CpxToWavelet2D( T ratio ) : m_ratio( ratio ) {}
  template <class Tuple>
  void operator()(Tuple in) const {
    T tmp2 = boost::get<2>(in);
    boost::get<2>(in) = (boost::get<1>(in)-boost::get<3>(in))*m_ratio;
    boost::get<1>(in) = (boost::get<1>(in)+boost::get<3>(in))*m_ratio;
    boost::get<3>(in) = (-boost::get<0>(in)+tmp2) * m_ratio;
    boost::get<0>(in) = (boost::get<0>(in)+tmp2) * m_ratio;
  }
  const T m_ratio;
};

/** \struct WaveletToCpx3D
 * \brief A linear operator that allows to perform the alternative y/z negative
 * frequency cancellation
 *
 * From the 8 raw output of the dtcwt filtering process, one actually obtains
 * a combination of signals(real): a , and their hilbert transform(imag): b
 * We recall that the complex resulting vector 1/2 (a + i b) has its negative
 * frequency content cancelled.
 * 1/2 (a - ib) has its positive frequency cancelled
 *
 * More precisely, let say we have 8 row vectors as input, they
 * would represent:
 * 0: real Z * real Y * real X : real
 * 1: real Z * real Y * Imag X : imag
 * 2: real Z * imag Y * real X : imag
 * 3: real Z * imag Y * Imag X : -real
 * 4: imag Z * real Y * real X : real
 * 5: imag Z * real Y * Imag X : imag
 * 6: imag Z * imag Y * real X : imag
 * 7: imag Z * imag Y * Imag X : -real
 *
 * When doing a combination of real and Hilbert/imaginary filters along x,y,z
 * where negative frequency cancellation property is used, we obtain a
 * combination that has a simple but limited interpretation:
 * (realZ+iImagZ)(realY+iImagY)(realX+iImagX)
 * =   realZrealYrealX - RealZImagYImagX+i(realZrealYImagX + realZImagYRealX)
 * + i(ImagZrealYrealX - ImagZImagYImagX)-(ImagZrealYImagX + ImagZImagYRealX)
 *
 * =  realZrealYrealX - RealZImagYImagX - ImagZrealYImagX - ImagZImagYRealX
 * +i(realZrealYImagX + realZImagYRealX + ImagZrealYrealX - ImagZImagYImagX)
 * =(0-3-5-6) + i(1+2+4-7) : positiveX, positiveY, positiveZ freq quadrant
 *
 * One can use a negative z cancellation, keeping x/y positive cancellation
 * (realZ-iImagZ)(realY+iImagY)(realX+iImagX):
 * =   realZrealYrealX - RealZImagYImagX +i (realZrealYImagX + realZImagYRealX)
 * - i(ImagZrealYrealX - ImagZImagYImagX) + (ImagZrealYImagX + ImagZImagYRealX)
 *
 * =  realZrealYrealX - RealZImagYImagX + ImagZrealYImagX + ImagZImagYRealX
 * +i(realZrealYImagX + realZImagYRealX - ImagZrealYrealX + ImagZImagYImagX)
 * =(0-3+5+6) + i(1+2-4+7) : positiveX, positiveY, positiveZ freq quadrant
 *
 * One can use a negative y cancellation, keeping x/z positive cancellation
 * (realZ+iImagZ)(realY-iImagY)(realX+iImagX)
 * =   realZrealYrealX + RealZImagYImagX+i(realZrealYImagX - realZImagYRealX)
 * + i(ImagZrealYrealX + ImagZImagYImagX)-(ImagZrealYImagX - ImagZImagYRealX)
 *
 * =  realZrealYrealX + RealZImagYImagX - ImagZrealYImagX + ImagZImagYRealX
 * +i(realZrealYImagX - realZImagYRealX + ImagZrealYrealX + ImagZImagYImagX)
 * =(0+3-5+6) + i(1-2+4+7) : positiveX, positiveY, positiveZ freq quadrant
 *
 * One can use a negative z/y cancellation, keeping x positive cancellation
 * (realZ-iImagZ)(realY-iImagY)(realX+iImagX)
 * =   realZrealYrealX + RealZImagYImagX+i(realZrealYImagX - realZImagYRealX)
 *  -i(ImagZrealYrealX + ImagZImagYImagX)+(ImagZrealYImagX - ImagZImagYRealX)
 *
 * =  realZrealYrealX + RealZImagYImagX + ImagZrealYImagX - ImagZImagYRealX
 * +i(realZrealYImagX - realZImagYRealX - ImagZrealYrealX - ImagZImagYImagX)
 * =(0+3+5-6) + i(1-2-4-7) : positiveX, positiveY, positiveZ freq quadrant
 *
 * Resulting operator can be written: <br>
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
 */
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
 * \brief Just the inverse of WaveletToCpx3D
 *
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
 */
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
