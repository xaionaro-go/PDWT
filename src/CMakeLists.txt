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
    #cuda/common.cu
    cuda/coeffContainer.cu
    cuda/constant.cu
    #cuda/haar.cu
    #cuda/nonseparable.cu
    #cuda/separable.cu
    #cuda/wt.cu
  )
  set(headers_cuda
  )
endif ()
################
#  /Cuda Part  #
################i

###############
#  Src files  #
###############
set(sources
  filters.cpp
  coeffContainer.cpp
  Wavelet.cpp
  Wavelet1D.cpp
  Wavelet2D.cpp
  Wavelet3D.cpp
  # cuda source files
  ${sources_cuda}
)

set(headers
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
  ${boost_deps})
set_property(TARGET ${OUTPUT_NAME}
  PROPERTY OUTPUT_NAME ${output_name})
install(TARGETS ${OUTPUT_NAME}
  LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY}
)
