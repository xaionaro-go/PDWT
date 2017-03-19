#ifndef COEFFCONTAINER_H
#define COEFFCONTAINER_H

//STL
#include <list>

/** \class CoeffContainer
 * \brief Generic interface that should be implemented by wavelet coefficient
 * containers of arbitrary dimension.
 * It has been designed to ease the access to pyramidal data decomposition
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class SubContainerT>
class CoeffContainer {
 public:
  /// Defaulted Constructor
  CoeffContainer()=default;

  /// Defaulted Destructor
  virtual ~CoeffContainer()=default;

  /**
   * Function that will take care of initializing the container, especially it
   * is in charge of allocatin memory
   *
   * \param[in] shape number of elements in each direction for the input
   * \param[in] level depth of the wavelet decomposition, 0 is no decomposition
   * \param[in] value value used to fill the input, 0 by default
   *
   * \return void 
   */
  virtual void Initialize(std::list<std::size_t> shape, std::size_t level,
      T value = 0);

protected:
  /// The pyramid containing Various stage of the DT
  SubContainerT m_coeff;
  
  /// The initial image shape
  std::list<std::size_t> m_shape;

  /// Scale depth of the wavelet decomposition
  std::size_t m_level;

  /// Shape of each level of coefficients
  std::list<std::list<std::size_t>> m_scaleShape;

  /// Size of each level of coefficients
  std::list<std::size_t> m_scaleSize;
 };

/** \class CoeffContainer1D
 * \brief Implementation of the CoeffContainer interface for the one
 * dimensional case
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class SubContainerT>
class CoeffContainer1D : public CoeffContainer<T,SubContainerT> {
 public:
  /// Defaulted Constructor
  CoeffContainer1D()=default;

  /// Defaulted Destructor
  virtual ~CoeffContainer1D()=default;

  /// Return the dimensionality of the container
  virtual std::size_t GetNbDimension() const override { return m_dimensions; };

 protected:
  static const std::size_t m_dimensions=1;
};

/** \class CoeffContainer2D
 * \brief Implementation of the CoeffContainer interface for the two
 * dimensional case
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class SubContainerT>
class CoeffContainer2D : public CoeffContainer<T,SubContainerT> {
 public:
  /// Defaulted Constructor
  CoeffContainer2D()=default;

  /// Defaulted Destructor
  virtual ~CoeffContainer2D()=default;

  /// Return the dimensionality of the container
  virtual std::size_t GetNbDimension() const override { return m_dimensions; };

 protected:
  static const std::size_t m_dimensions=2;
};

/** \class CoeffContainer3D
 * \brief Implementation of the CoeffContainer interface for the three
 * dimensional case
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class SubContainerT>
class CoeffContainer3D : public CoeffContainer<T,SubContainerT> {
 public:
  /// Defaulted Constructor
  CoeffContainer3D()=default;

  /// Defaulted Destructor
  virtual ~CoeffContainer3D()=default;
 
  /// Return the dimensionality of the container 
  virtual std::size_t GetNbDimension() const override { return m_dimensions; };

 protected:
  static const std::size_t m_dimensions=3;
};

#endif /* COEFFCONTAINER_H */
