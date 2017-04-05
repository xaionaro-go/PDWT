#ifndef COEFFCONTAINER_H
#define COEFFCONTAINER_H

// STL
#include <algorithm>
#include <functional>
#include <vector>

// Local


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

  /// Allocating constructor
  CoeffContainer(std::vector<size_t> size, int nlevel) {
    Initialize(size, nlevel);
  }

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
  virtual void Initialize(std::vector<std::size_t> shape, std::size_t nlevel,
      T value = 0) {
    m_level = nlevel;
    std::vector<std::size_t> curShape{shape};
    auto sb = curShape.begin();
    auto se = curShape.end();
    auto divider = [](size_t in) {return in/2;};

    for (size_t i=0; i<=m_level; i++) {
      m_scaleSize.emplace_back(std::accumulate(
        sb,se,1u, std::multiplies<size_t>()));
      m_scaleShape.emplace_back(curShape);
      std::transform(sb,se,sb,divider);
    }
    // Allocate memory
    m_coeff.resize(m_scaleSize.front());
  }

  /// Return the dimensionality of the container
  virtual std::size_t GetNbDimension() const = 0;

  /// Returns a pointer to the Low frequency subspace for the given scale
  T* GetLowSubspacePtr(size_t scale, size_t band=0) {
    // Band offset, accounts for redundancy
	size_t bandOffset = band*m_scaleSize.at(0);
	// For each band, we have (2^d)-1 high frequency subband + 1 low frequency
    // We decided to put the low frequency at the end in our layout
	size_t nbHighFreqSubband = std::pow(2,GetNbDimension())-1;
	size_t scaleOffset = nbHighFreqSubband*std::accumulate(
      m_scaleSize.begin()+1, m_scaleSize.begin()+scale+2,
      0u, std::plus<size_t>());
	return m_coeff.data()+bandOffset+scaleOffset;
  }

  /// Return a pointer to the High frequency subspace for the given level
  T* GetHighSubspacePtr(size_t scale, size_t subband, size_t band=0) {
	// Band offset, accounts for redundancy
	size_t bandOffset = band*m_scaleSize.at(0);
	// For each band, we have (2^d)-1 high frequency subband + 1 low frequency
    // We decided to put the low frequency at the end in our layout
	size_t nbHighFreqSubband = std::pow(2,GetNbDimension())-1;
	size_t scaleOffset = nbHighFreqSubband*std::accumulate(
      m_scaleSize.begin()+1, m_scaleSize.begin()+scale+1,
      0u, std::plus<size_t>());
	// subband index accounts for filtering combination, ie ranges from 0 to
    // (2^d)-2
	size_t subbandOffset = subband*m_scaleSize.at(scale+1);
	return  m_coeff.data()+bandOffset+scaleOffset+subbandOffset;
  }

 protected:
  /// The pyramid containing Various stage of the DT
  SubContainerT m_coeff;
  
  /// Scale depth of the wavelet decomposition
  std::size_t m_level;

  /// Shape of each level of coefficients
  std::vector<std::vector<std::size_t>> m_scaleShape;

  /// Size of each level of coefficients
  std::vector<std::size_t> m_scaleSize;
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

  /// Allocating constructor
  CoeffContainer1D(std::vector<size_t> size, int nlevel) :
    CoeffContainer<T,SubContainerT>(size,nlevel) {}

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

/** \class CoeffContainerCpx
 * \brief Implementation of the CoeffContainer interface for storing complex
 * wavelet transform, for the generic case (arbitrary dimension)
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class SubContainerT>
class CoeffContainerCpx : public CoeffContainer<T,SubContainerT> {
 public:
  /// Deleted Constructor
  CoeffContainerCpx()=delete;

  /// Defaulted Destructor
  virtual ~CoeffContainerCpx()=default;
};

#endif /* COEFFCONTAINER_H */
