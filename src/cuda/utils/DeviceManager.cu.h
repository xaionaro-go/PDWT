#ifndef 
#define DEVICEMANAGER_CU_H_

// STL
#include <vector>

// Cuda
#include <cuda_runtime.h>

#define CUDA_MINIMUM_COMPUTE_CAPABILITY_MAJOR 3

/** \struct deviceDesc
 * \brief gives a short description of each device, along with accessors.
 *
 * \author Thibault Notargiacomo
 */
struct DeviceDesc {
  int id; //TODO TN check if id is not already contained in CudaDeviceProp
  cudaDeviceProp deviceProp;
  std::vector<cudaStream_t> vStream;
};


/** \class deviceManager
 * \brief List and allow access to multiple gpu device on the local computer.
 *
 * This singleton list all available GPUs at startup of the application.
 * Then it allows to access each gpu.
 *
 * \author Thibault Notargiacomo
 */
class DeviceManager {
 public:
  static DeviceManager& GetInstance();

  void reset(size_t nbDevice = 0);
  const std::vector<DeviceDesc>& GetDeviceDesc() const;
  int GetCurrentDeviceIdx() const;

 protected:
  virtual ~DeviceManager();
  DeviceManager();
  std::vector<DeviceDesc> m_vDeviceDesc;
};

#endif // DEVICEMANAGER_CU_H_
