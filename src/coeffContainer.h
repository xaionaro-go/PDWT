#ifndef COEFFCONTAINER_H
#define COEFFCONTAINER_H

// STL
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <memory>
#include <numeric>
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
  CoeffContainer(const std::vector<size_t>& size, int nlevel) {
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
  virtual void Initialize(const std::vector<size_t>& shape, size_t nlevel,
      T value = 0) {
    m_level = nlevel;
    std::vector<size_t> curShape{shape};
    auto sb = curShape.begin();
    auto se = curShape.end();
    auto divider = [](size_t in) {return (in+(in&1))/2;};

    for (size_t i=0; i<=m_level; i++) {
      size_t levelSize = std::accumulate(sb,se,1u, std::multiplies<size_t>());
      assert(levelSize>0);
      m_scaleSize.emplace_back(levelSize);
      m_scaleShape.emplace_back(curShape);
      std::transform(sb,se,sb,divider);
    }
    /** Compute total size in the coeff space.
     * It is equal to the sum of the size of all detail subspaces,
     * plus the size of the last approximation space
     */
    size_t totalSize=std::accumulate(m_scaleSize.cbegin()+1,
      m_scaleSize.cend(),0u)+m_scaleSize.back();
    // Allocate memory
    m_coeff.resize(totalSize);
    // Allocate temporary buffer
    if (m_level>1) {
      m_vptcoeff.push_back(std::make_unique<SubContainerT>(m_scaleSize.at(1)));
      m_vptcoeff.push_back(std::make_unique<SubContainerT>(m_scaleSize.at(2)));
    }
  }

  /// Return the dimensionality of the container
  virtual size_t GetNbDimension() const = 0;

  /// Returns a pointer to the Low frequency subspace for the given scale
  T* GetLowSubspacePtr(size_t scale, size_t band=0) {
    // Band offset, accounts for redundancy
	size_t bandOffset = band*m_scaleSize.at(0);
	// For each band, we have (2^d)-1 high frequency subband + 1 low frequency
    // We decided to put the low frequency at the end in our layout
	size_t nbHighFreqSubband = std::pow(2,GetNbDimension())-1;
    assert(m_scaleSize.size()>1);
	size_t scaleOffset = nbHighFreqSubband*std::accumulate(
      m_scaleSize.cbegin()+1,
      std::min(m_scaleSize.cbegin()+scale+2,m_scaleSize.cend()),
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

  /// Return a pointer to a temporary buffer
  std::vector<std::unique_ptr<SubContainerT>>& GetTmpBuffPtr() {
    return m_vptcoeff;
  }

  /// Simple proxy for subcontainer begin iterator getter
  auto begin() { return m_coeff.begin(); }
  /// Simple proxy for subcontainer end iterator getter
  auto end() { return m_coeff.end(); }
  /// Simple proxy for subcontainer const begin iterator getter
  auto cbegin() const { return m_coeff.cbegin(); }
  /// Simple proxy for subcontainer const end iterator getter
  auto cend() const { return m_coeff.cend(); }

  /// Returns the size of data for the given scale
  auto GetScaleShape(size_t scale) {
    return m_scaleShape.at(scale);
  }

 protected:
  /// The pyramid containing Various stage of the DT
  SubContainerT m_coeff;

  /// A temporary buffer used to compute DWT
  std::vector<std::unique_ptr<SubContainerT>> m_vptcoeff;  

  /// Scale depth of the wavelet decomposition
  size_t m_level;

  /// Shape of each level of coefficients
  std::vector<std::vector<size_t>> m_scaleShape;

  /// Size of each level of coefficients
  std::vector<size_t> m_scaleSize;
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
  virtual size_t GetNbDimension() const override { return m_dimensions; };

 protected:
  static const size_t m_dimensions=1;
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
  virtual size_t GetNbDimension() const override { return m_dimensions; };

 protected:
  static const size_t m_dimensions=2;
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
  virtual size_t GetNbDimension() const override { return m_dimensions; };

 protected:
  static const size_t m_dimensions=3;
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
