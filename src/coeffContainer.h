#ifndef COEFFCONTAINER_H
#define COEFFCONTAINER_H

// STL
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

// Boost
#include <boost/iterator/zip_iterator.hpp>

// Local
#include "dualTreeLinearOp.h"

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
  virtual void InitializeSizes(const std::vector<size_t>& shape, size_t nlevel,
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
     * times the number of band if one uses dual tree scheme
     */
    m_nbSubBand = std::pow(2,GetNbDimension())-1;
    m_nbBand = DoUseDualTreeBand()?std::pow(2,GetNbDimension()):1;
    m_bandSize = m_nbSubBand * std::accumulate(m_scaleSize.cbegin()+1,
      m_scaleSize.cend(),0u) + m_scaleSize.back();
    assert(m_scaleSize.size()>0);

    //std::cout <<"nbBand is "<<nbBand<<" nbSubBand is "<<nbSubBand<<std::endl;
    //auto print = [](auto in) {std::cout<<" , " <<in;};
    //std::cout<<"scaleSize is: ";
    //std::for_each(m_scaleSize.cbegin(),m_scaleSize.cend(),print);
    //std::cout<<std::endl;
    //std::cout<<"Total size is "<<totalSize<<"x"<<std::endl;
  }

  /// Return the dimensionality of the container
  virtual size_t GetNbDimension() const = 0;
  
  /// Return wether we are using the dual tree scheme or not
  virtual bool DoUseDualTreeBand() const { return false; };

  /// Returns a pointer to the Low frequency subspace for the given scale
  T* GetLowSubspacePtr(size_t scale, size_t band=0) {
    // Band offset, accounts for redundancy
	size_t bandOffset = band*m_bandSize;
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
	size_t bandOffset = band*m_bandSize;
	// For each band, we have (2^d)-1 high frequency subband + 1 low frequency
    // We decided to put the low frequency at the end in our layout
	size_t nbHighFreqSubband = std::pow(2,GetNbDimension())-1;
	size_t scaleOffset = nbHighFreqSubband*std::accumulate(
      m_scaleSize.begin()+1, m_scaleSize.begin()+scale+1,
      0u, std::plus<size_t>());
	// subband index accounts for filtering combination, ie ranges from 0 to
    // (2^d)-2
	size_t subbandOffset = subband*m_scaleSize.at(scale+1);
    //std::cout<<"Band offset is "<<bandOffset<<std::endl;
    //std::cout<<"scaleOffset is "<<scaleOffset<<std::endl;
    //std::cout<<"subband offset is "<<subbandOffset<<std::endl;
	return  m_coeff.data()+bandOffset+scaleOffset+subbandOffset;
  }

  /// Apply functor to all complex coefficients (as a tuple)
  template<class Func>
  void ApplyCpxFunctor(Func cpxFunctor) {
    auto dtBegin = m_coeff.begin();

    for(size_t bandIdx=0; bandIdx<m_nbBand; bandIdx++) {

      auto bandBegin = dtBegin+2*bandIdx*m_bandSize;
      //Build a real/imag pair-iterator using zip
      auto cpxBegin(boost::make_zip_iterator(boost::make_tuple(
        bandBegin, bandBegin+m_bandSize)));
      auto cpxEnd(boost::make_zip_iterator(boost::make_tuple(
        bandBegin+m_bandSize, bandBegin+2*m_bandSize)));

      std::for_each(cpxBegin, cpxEnd, cpxFunctor);
    }
  }

  /// Return a pointer to a temporary buffer
  virtual T* GetTmpBuffPtr(size_t level, size_t index)=0;

  /// Simple proxy for subcontainer begin iterator getter
  auto begin() { return m_coeff.begin(); }
  /// Simple proxy for subcontainer end iterator getter
  auto end() { return m_coeff.end(); }
  /// Simple proxy for subcontainer const begin iterator getter
  auto cbegin() const { return m_coeff.cbegin(); }
  /// Simple proxy for subcontainer const end iterator getter
  auto cend() const { return m_coeff.cend(); }
  /// Get size of the underlying container
  size_t size() const { return m_coeff.size(); };

  /// Returns the size of data for the given scale
  auto GetScaleShape(size_t scale) {
    return m_scaleShape.at(scale);
  }

 protected:
  // Allocate main memory
  virtual void AllocateMainBuffer() {
    size_t totalSize = m_nbBand * m_bandSize;
    m_coeff.resize(totalSize);
  }

  // Allocate temporary buffer
  virtual void AllocateTmpBuffer()=0;

 protected:
  /// The pyramid containing Various stage of the DT
  SubContainerT m_coeff;

  /// A temporary buffer used to compute DWT
  std::unique_ptr<SubContainerT> m_ptcoeff;  

  /// Scale depth of the wavelet decomposition
  size_t m_level;

  /// Shape of each level of coefficients
  std::vector<std::vector<size_t>> m_scaleShape;

  /// Size of each level of coefficients
  std::vector<size_t> m_scaleSize;

  /**
  * Number of subband: only depends on the dimension of the data for
  * separable wavelets
  **/
  size_t m_nbSubBand;
  
  /// Number of band in case of dual/quad/oct-tree wavelet
  size_t m_nbBand;
  
  /// Size of each band in case of dual/quad/oct-tree wavelet
  size_t m_bandSize; 
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
  CoeffContainer1D(std::vector<size_t> size, int nlevel) {
    this->InitializeSizes(size, nlevel);
    // Allocate memory
    this->AllocateMainBuffer();
    // Allocate temporary buffer
    this->AllocateTmpBuffer();
  }

  /// Defaulted Destructor
  virtual ~CoeffContainer1D()=default;

  /// Return the dimensionality of the container
  virtual size_t GetNbDimension() const override { return m_dimensions; };

  /// Return a pointer to a temporary buffer
  virtual T* GetTmpBuffPtr(size_t level, size_t bandIdx=0) override {
    if (this->m_level>1) {
      // Layout is repeated nbBand folds, and contains every time
      // scaSize[1]+scaleSize[2]
      size_t bandOffset = this->m_scaleSize.at(1)+this->m_scaleSize.at(2);
      /**
       * Temp buffers are designed such that there are successive swap, and
       * we have the following behaviour:
       * For forward transform:
       * When starting at level 0, we need a temp buffer of size scalesize[1]
       * For backward transform:
       * when arriving at level 2, we can put the reconstruction of level 1 in
       * a temporary buffer of size scalesize[1]
       */
      size_t levelOffset = (level%2)*this->m_scaleSize.at(1);
      return this->m_ptcoeff->data()+bandOffset*bandIdx+levelOffset;
    } else {
      assert(false);
      return nullptr;
    }
  }

 protected:
  // Allocate temporary buffer
  virtual void AllocateTmpBuffer() override {
    if (this->m_level>1) {
      this->m_ptcoeff=std::make_unique<SubContainerT>(
        this->m_nbBand*(this->m_scaleSize.at(1)+this->m_scaleSize.at(2)));
    }
  }
 protected:
  static const size_t m_dimensions=1;
};

/** \class DTCoeffContainer1D
 * \brief Implementation of the CoeffContainer interface for the one
 * dimensional case, compatible with the dual tree scheme
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class SubContainerT>
class DTCoeffContainer1D : public CoeffContainer1D<T,SubContainerT> {
 public:
  /// Defaulted Constructor
  DTCoeffContainer1D()=default;

  /// Allocating constructor
  DTCoeffContainer1D(std::vector<size_t> size, int nlevel) {
    this->InitializeSizes(size, nlevel);
    // Allocate memory
    this->AllocateMainBuffer();
    // Allocate temporary buffer
    this->AllocateTmpBuffer();
  }

  /// Defaulted Destructor
  virtual ~DTCoeffContainer1D()=default;

  /// Return wether we are using the dual tree scheme or not
  virtual bool DoUseDualTreeBand() const override { return true; };

  /// Make the magical mixture of negative frequency cancelling signals
  int WaveletToCpx() {
    std::transform(this->m_coeff.begin(),this->m_coeff.end(),
      this->m_coeff.begin(),
      [](auto in) { return in*m_normalizationRatio; });
    return 1;
  }

  /// Should be the exact inverse mapping of WaveletToCpx
  int CpxToWavelet() { return WaveletToCpx(); };

 protected:
  static constexpr const T m_normalizationRatio = 1.0/std::sqrt(2);
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

  /// Allocating constructor
  CoeffContainer2D(std::vector<size_t> size, int nlevel) {
    Initialize2DContainer(size, nlevel);
  }

  /// Defaulted Destructor
  virtual ~CoeffContainer2D()=default;

  /// Initializes both sizes and buffers
  void Initialize2DContainer(std::vector<size_t>& size, int nlevel) {
    this->InitializeSizes(size, nlevel);
    /**
     * tmp buffer should be able to store output of single Y filtering,
     * size is initiale size where only second dimension has been divided
     * by 2
     */
    m_tmpSingleBuffSize=this->m_scaleShape.at(0).at(0)*
      this->m_scaleShape.at(1).at(1);
    if (this->m_level>1) {
      m_tmpBuffBandOffset=m_tmpSingleBuffSize*2+this->m_scaleSize.at(1);
    } else {
      m_tmpBuffBandOffset=m_tmpSingleBuffSize*2;
    }
    // Allocate memory
    this->AllocateMainBuffer();
    // Allocate temporary buffer
    this->AllocateTmpBuffer();
  }

  /// Return the dimensionality of the container
  virtual size_t GetNbDimension() const override { return m_dimensions; };

  /// Return a pointer to a temporary buffer, idx stands for low(0) or high(1)
  virtual T* GetHalfTmpBuffPtr(size_t subBandIdx, size_t bandIdx=0) {
    // Layout is  described in constructor
    // subbandoffset
    size_t subBandOffset = m_tmpSingleBuffSize;
    return this->m_ptcoeff->data()+m_tmpBuffBandOffset*bandIdx+
      subBandOffset*subBandIdx;
  }

  /**
   * Return a pointer to a temporary buffer, for lowpass tmp storage
   * The rationale behind this choice is the following:
   * For forward transform:
   * Y filtering: Worker has to write lowpass at first slot and highpass at
   * second slot
   * X filtering: Worker has to write lowpass at OutLowTmpBuffer and Highpass
   * to wavelt tree, while reading both slot1 and slot2 fron tmp
   * For backward transform:
   * We have two independant X filtering to invert:
   * LowYLowX and LowYHighX goes into first slot
   * HighYLowX and HighYHighX goes into second slot
   * then reconstruction can go into the third LowTmpBuff, or image level>=1 
   */
  virtual T* GetOutLowTmpBuffPtr(size_t bandIdx=0) {
    // Layout is defined in constructor
    size_t lowPassOffset = 2*m_tmpSingleBuffSize;
    return this->m_ptcoeff->data()+m_tmpBuffBandOffset*bandIdx+lowPassOffset;
  }

  /// Return a pointer to a temporary buffer
  //TODO TN: a refaire
  virtual T* GetTmpBuffPtr(size_t level, size_t index) {
    return GetHalfTmpBuffPtr(0,index);
  }

 protected:
  // Allocate temporary buffer
  virtual void AllocateTmpBuffer() override {
    this->m_ptcoeff=std::make_unique<SubContainerT>(
      this->m_nbBand*m_tmpBuffBandOffset);
  }

 protected:
  size_t m_tmpSingleBuffSize;
  size_t m_tmpBuffBandOffset;
  static const size_t m_dimensions=2;
};

/** \class DTCoeffContainer2D
 * \brief Implementation of the CoeffContainer interface for the two
 * dimensional case, compatible with the dual tree scheme
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class SubContainerT>
class DTCoeffContainer2D : public CoeffContainer2D<T,SubContainerT> {
 public:
  /// Defaulted Constructor
  DTCoeffContainer2D()=default;

  /// True constructor
  DTCoeffContainer2D(std::vector<size_t> size, int nlevel) {
    this->Initialize2DContainer(size, nlevel);
  }

  /// Defaulted Destructor
  virtual ~DTCoeffContainer2D()=default;

  /// Return wether we are using the dual tree scheme or not
  virtual bool DoUseDualTreeBand() const override { return true; };

  /// Apply functor to all complex coefficients (as a tuple)
  template<class Func>
  void ApplyBandFunctor(const Func& bFunctor) {
    auto dtBegin = this->m_coeff.begin();
    size_t bSize = this->m_bandSize;

    //Build a band 0/1/2/3 4-uplet iterator using zip
    auto bBegin = boost::make_zip_iterator(boost::make_tuple(
      dtBegin,
      dtBegin+bSize,
      dtBegin+2*bSize,
      dtBegin+3*bSize));
    auto bEnd = boost::make_zip_iterator(boost::make_tuple(
      dtBegin+bSize,
      dtBegin+2*bSize,
      dtBegin+3*bSize,
      dtBegin+4*bSize));

    std::for_each(bBegin, bEnd, bFunctor);
  }

  /// Make the magical mixture of negative frequency cancelling signals
  int WaveletToCpx() {
    ApplyBandFunctor(WaveletToCpx2D<T>(m_normalizationRatio));
    return 1;
  }

  /// Should be the exact inverse mapping of WaveletToCpx
  int CpxToWavelet() {
    ApplyBandFunctor(CpxToWavelet2D<T>(m_normalizationRatio));
    return 1;
  }

 protected:
  static constexpr const T m_normalizationRatio = 1.0/(2.0*std::sqrt(2));
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

  // True constructor
  CoeffContainer3D(std::vector<size_t> size, int nlevel) {
    Initialize3DContainer(size, nlevel);
  }

  /// Defaulted Destructor
  virtual ~CoeffContainer3D()=default;
 
  /// Initializes both sizes and buffers
  void Initialize3DContainer(std::vector<size_t>& size, int nlevel) {
    this->InitializeSizes(size, nlevel);
    /**
     * in 3D wavelet transform, memory consumption is critical, this is why
     * we will perform all steps sequentially.
     *
     * FORWARD
     * Step 1:
     * tmp buffer should be able to store output of one single Z filtering,
     * which is initial size where only third dimension has been divided
     * by 2: SizeX*SizeY*(SizeZ/2)
     * Step 2:
     * We then need to perform Y filtering, and store one single output of
     * size initial size with dimension 2 and 3 divided by 2:
     * SizeX*(SizeY/2)*(SizeZ/2)
     * Step 3
     * For X filtering, data can be written directly to coefficient storage
     * but we also need space for the lowpass output for next transform if
     * there is more than one level. In addition input/output need to be non
     * overlapping, so that adds an additional space for next output if needed
     * worst case needed size is
     * (SizeX/2)*(SizeY/2)*(SizeZ/2)+(SizeX/4)*(SizeY/4)*(SizeZ/4)
     *
     * BACKWARD
     * Step 1, multiple sequential pass: foreach Ysubsp, foreach Z subsp:
     * tmp buffer stores the output of a double (l/h) x filtering of size
     * (worst case:) SizeX*(SizeY/2)*(SizeZ/2)
     * Step 2:
     * The xfiltered data is then filtered along Y, with a single filter, and
     * written (or in an update manner if not first) to a second tmp buffer
     * (worst case:) SizeX*SizeY*(SizeZ/2)
     * Step 3:
     * Either we are at the last reconstruction step, and one doesn't need
     * any additional memory (writing in output directly) or the worst case
     * then need to write the low pass approximation of size 
     * (SizeX/2)*(SizeY/2)*(SizeZ/2) somewhere, from another buffer of size
     * (SizeX/4)*(SizeY/4)*(SizeZ/4)
     * 
     * One can conclude that Forward transform needs more tmp memory
     *
     * TODO TN: DT extension (very bad)
     * Z filter and Y filter mandatory buff are always allocated in one single
     * instance. However, Z tmp output and Z_2 tmp output are allocated
     * in nbBand instances
     */
    m_tmpZOutSingleSize=this->m_scaleShape.at(0).at(0)*
      this->m_scaleShape.at(0).at(1)*
      this->m_scaleShape.at(1).at(2);
    m_tmpYOutSingleSize=this->m_scaleShape.at(0).at(0)*
      this->m_scaleShape.at(1).at(1)*
      this->m_scaleShape.at(1).at(2);
    m_tmpXOutSingleSize=this->m_scaleSize.at(1);

    m_tmpBuffBaseOffset=m_tmpZOutSingleSize+m_tmpYOutSingleSize;
    if (this->m_level>2) {
      m_tmpBuffBandOffset=m_tmpXOutSingleSize+this->m_scaleSize.at(2);
    } else if (this->m_level>1) {
      m_tmpBuffBandOffset=m_tmpXOutSingleSize;
    } else {
      m_tmpBuffBandOffset=0;
    }
    // Allocate memory
    this->AllocateMainBuffer();
    // Allocate temporary buffer
    this->AllocateTmpBuffer();
  }

  /// Return the dimensionality of the container
  virtual size_t GetNbDimension() const override { return m_dimensions; };

  /// Return a pointer to a temporary buffer, idx stands for low(0) or high(1)
  virtual T* GetHalfTmpBuffPtr(size_t subBandIdx, size_t bandIdx) {
    // Layout is  described in constructor
    // TODO TN: for TN, layout is very bad, and should be designed properly
    size_t subBandOffset=0;
    size_t fullBuffTmp=0;
    if (subBandIdx==1) {
      subBandOffset = m_tmpZOutSingleSize;
    } else if(subBandIdx==2) {
      subBandOffset = m_tmpBuffBaseOffset;
    } else if(subBandIdx==3) {
      subBandOffset = m_tmpBuffBaseOffset;
      fullBuffTmp = m_tmpXOutSingleSize;
    }
    return this->m_ptcoeff->data()+subBandOffset+bandIdx*m_tmpBuffBandOffset+
      fullBuffTmp;
  }

  /**
   * Return a pointer to a temporary buffer, for lowpass tmp storage
   */
  virtual T* GetOutLowTmpBuffPtr(size_t bandIdx=0) {
    // Layout is defined in constructor
    return GetHalfTmpBuffPtr(2, bandIdx);
  }

  /// Return a pointer to a temporary buffer
  //TODO TN: a refaire
  virtual T* GetTmpBuffPtr(size_t level, size_t index) {
    return nullptr;
  }

 protected:
  // Allocate temporary buffer
  virtual void AllocateTmpBuffer() override {
    this->m_ptcoeff=std::make_unique<SubContainerT>(
      m_tmpBuffBaseOffset+this->m_nbBand*m_tmpBuffBandOffset);
  }

 protected:
  size_t m_tmpBuffBaseOffset;
  size_t m_tmpBuffBandOffset;
  size_t m_tmpZOutSingleSize;
  size_t m_tmpYOutSingleSize;
  size_t m_tmpXOutSingleSize;
  static const size_t m_dimensions=3;
};

/** \class DTCoeffContainer3D
 * \brief Implementation of the CoeffContainer interface for the three
 * dimensional case, compatible with the dual tree scheme
 *
 * \author Thibault Notargiacomo
 */
template<typename T, class SubContainerT>
class DTCoeffContainer3D : public CoeffContainer3D<T,SubContainerT> {
 public:
  /// Defaulted Constructor
  DTCoeffContainer3D()=default;

  /// Allocating constructor
  DTCoeffContainer3D(std::vector<size_t> size, int nlevel) {
    this->InitializeSizes(size, nlevel);
    // Allocate memory
    this->AllocateMainBuffer();
    // Allocate temporary buffer
    this->AllocateTmpBuffer();
  }

  /// Defaulted Destructor
  virtual ~DTCoeffContainer3D()=default;

  /// Return wether we are using the dual tree scheme or not
  virtual bool DoUseDualTreeBand() const override { return true; };

  /// Apply functor to all complex coefficients (as a tuple)
  template<class Func>
  void ApplyBandFunctor(const Func& bFunctor) {
    auto dtBegin = this->m_coeff.begin();
    size_t bSize = this->m_bandSize;

    //Build a band 0/1/2/3 4-uplet iterator using zip
    auto bBegin = boost::make_zip_iterator(boost::make_tuple(
      dtBegin,
      dtBegin+bSize,
      dtBegin+2*bSize,
      dtBegin+3*bSize,
      dtBegin+4*bSize,
      dtBegin+5*bSize,
      dtBegin+6*bSize,
      dtBegin+7*bSize));
    auto bEnd = boost::make_zip_iterator(boost::make_tuple(
      dtBegin+bSize,
      dtBegin+2*bSize,
      dtBegin+3*bSize,
      dtBegin+4*bSize,
      dtBegin+5*bSize,
      dtBegin+6*bSize,
      dtBegin+7*bSize,
      dtBegin+8*bSize));

    std::for_each(bBegin, bEnd, bFunctor);
  }

  /// Make the magical mixture of negative frequency cancelling signals
  int WaveletToCpx() {
    ApplyBandFunctor(WaveletToCpx3D<T>(m_normalizationRatio));
    return 1;
  }

  /// Should be the exact inverse mapping of WaveletToCpx
  int CpxToWavelet() {
    ApplyBandFunctor(CpxToWavelet3D<T>(m_normalizationRatio));
    return 1;
  }

 protected:
  static constexpr const T m_normalizationRatio = 1.0/(4.0*std::sqrt(2));
};

#endif /* COEFFCONTAINER_H */
