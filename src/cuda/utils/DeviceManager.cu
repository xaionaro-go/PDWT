//Local
#include "DeviceManager.cu.h"

//STL
#include <algorithm>
#include <cassert>
#include <numeric>
#include <iostream>
#include <iomanip>

//OMP
#include <omp.h>

#define NB_BYTE_PER_GIGA 1024*1024*1024

DeviceManager& CudaMultiGPUManager::GetInstance() {
  static DeviceManager instance;
  return instance;
}

void DeviceManager::reset(size_t nbDevice) {
  //clear current device list
  m_vDeviceDesc.clear();

  //Use cuda API to count number of devices
  int maxNbDevice = 0;
  checkCudaErrors(cudaGetDeviceCount(&maxNbDevice));

  //Build device list
  for (int devId=0; devId<maxNbDevice; devId++) {
    checkCudaErrors(cudaSetDevice(devId));
    checkCudaErrors(cudaDeviceReset());

    checkCudaErrors(cudaSetDeviceFlags(cudaDeviceBlockingSync));

    cudaDeviceProp devProp;
    bool bComputeCompatibility=true;
    checkCudaErrors(cudaGetDeviceProperties(&devProp,devId));

    if (devProp.major<=MIN_CUDA_COMPUTE_CAPABILITY_MAJOR) {
      std::cout<<"Device "<<devId<<" compute capability is too old"<<std::endl;
      //TODO TN: cudaDeviceGetByPCIBusId()
      bComputeCompatibility=false;
    }

    if( (m_vDeviceDesc.size()<nbDevice || nbDevice<=0) &&
        bComputeCompatibility) {
      checkCudaErrors(cudaSetDevice(devId));
      //Build DeviceDesc
      DeviceDesc devDesc;
      devDesc.id=devId;
      devDesc.deviceProp=devProp
      //By default, create only one stream
      cudaStream_t cudaStream;
      checkCudaErrors(cudaStreamCreate(&cudaStream));
      devDesc.vStream.push_back(cudaStream);
      //Finally push back new device descriptor in vector
      m_vDeviceDesc.push_back(devDesc);
    }
  }
}

const std::vector<DeviceDesc>& DeviceManager::GetDeviceDesc() const {
  return m_vDeviceDesc;
}

int DeviceManager::GetCurrentDeviceIdx() const {
  int devId;
  checkCudaErrors(cudaGetDevice(&devId));

  auto it=std::find(m_vDeviceDesc.cbegin(), m_vDeviceDesc.cend(),
    [](auto in) {return in.id==devId;} );

  assert(it!=m_vDeviceDesc.cend() );
  return std::distance(m_vDeviceDesc.cbegin(), it) ;
}

DeviceManager::~DeviceManager() {
  for (auto& it : m_vDeviceDesc ) {
    for (auto& stream : it.vStream) {
      checkCudaErrors(cudaStreamDestroy(stream));
    }
  }
}

DeviceManager::DeviceManager() {
  reset();

  std::cout << "******************************"
    "**********************************"<<std::endl;
  std::cout << "*  ID                          "
    "Name                   Memory   *"<<std::endl;
  for (const auto& dev : m_vDeviceDesc) {
    std::cout<<"*  "<<dev.id<<std::setw(31)<<
      dev.deviceProp.name<<std::setw(25)<<
      dev.deviceProp.totalGlobalMem/NB_BYTE_PER_GIGA<<"Go *"<<std::endl;
  }
  std::cout << "*                              "
    "    *"<<std::endl;
  std::cout << "*******************************"
    "*********************************"<<std::endl;
}
