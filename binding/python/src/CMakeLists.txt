cmake_minimum_required (VERSION 3.5)

##################
#  Generic Part  #
##################
project(pyPDWT)
set(OUTPUT_NAME pyPDWT)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Compiler must support some option for pybind11 code optimization
if (UNIX)
  add_definitions(${CMAKE_CXX_FLAGS} "-flto -fvisibility=hidden -std=c++14")
endif ()

#########################
#  Python binding code  #
#########################
set(sources
  CoeffContainer1D.cpp
  Wavelet1D.cpp
#  Wavelet2D.cpp
#  Wavelet3D.cpp
)

set(headers
)

####################
#  Lib Dependancy  #
####################
include_directories("${PRJ_ROOT_DIR}/src")

##############################
#  Equivalent of add_lirary  #
##############################
pybind11_add_module(${OUTPUT_NAME} ${sources})
