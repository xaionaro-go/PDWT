cmake_minimum_required(VERSION 3.5)

#################
#  Generic Part #
#################
project(pdwt)
set(OUTPUT_NAME pdwt)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

###############
#  Cuda Part  #
###############
if (CUDA_FOUND AND USE_CUDA)
  set(sources_cuda
    cuda/CoeffContainer.cu
    cuda/Constant.cu

    #utils
    cuda/utils/DeviceManager.cu
  )
  set(headers_cuda
    cuda/DeviceManager.cu.h
    cuda/ManagedAllocator.cu.h

    #utils
    cuda/utils/DeviceManager.cu.h
  )
endif ()
################
#  /Cuda Part  #
################i

###############
#  Src files  #
###############
set(sources
  Filters.cpp
  CoeffContainer.cpp
  Wavelet.cpp
  Wavelet1D.cpp
  Wavelet2D.cpp
  Wavelet3D.cpp
  # cuda source files
  ${sources_cuda}
)

set(headers
  Filters.h
  Wavelet.h
  Wavelet1D.h
  Wavelet2D.h
  Wavelet3D.h
  # cuda header files
  ${headers_cuda}
)

################
#  Executable  #
################

if (CUDA_FOUND)
  cuda_add_library(${OUTPUT_NAME} SHARED ${sources})
  CUDA_ADD_CUBLAS_TO_TARGET(${OUTPUT_NAME})
else ()
  add_library(${OUTPUT_NAME} SHARED ${sources})
endif ()
target_link_libraries(${OUTPUT_NAME}
  ${cuda_deps}
  ${mpi_deps}
  ${boost_deps}
  ${openmp_deps})
set_property(TARGET ${OUTPUT_NAME}
  PROPERTY OUTPUT_NAME ${output_name})
install(TARGETS ${OUTPUT_NAME}
  LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY}
)
