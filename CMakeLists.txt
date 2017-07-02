cmake_minimum_required(VERSION 3.5)

project(PDWT)
set(PROJECT_NAME PDWT)

#Default build destination
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(DOC_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc")
#Explicit make install destination
if (INSTALLDIR)
  set(RUNTIME_OUTPUT_DIRECTORY ${INSTALL_DIR}/bin)
  set(LIBRARY_OUTPUT_DIRECTORY ${INSTALL_DIR}/lib)
else ()
  set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
  set(LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
endif ()

#cmake ref
set(PRJ_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${PRJ_CMAKE_DIR})

#local ref
set(PRJ_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#################################
#  Generic Compilation options  #
#################################

#Compiler must support c++14
if (UNIX)
  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++14")
  add_definitions(${CMAKE_CXX_FLAGS_DEBUG} "-Wall -O0 -ggdb")
  add_definitions(${CMAKE_CXX_FLAGS_RELEASE} "-O3 -mtune=native -march=native")
endif ()

##############################
#  Cuda Compilation options  #
##############################

find_package(CUDA 8.0)

if (CUDA_FOUND AND USE_CUDA)
  add_definitions(-DUSE_CUDA)
  
  if (THRUST_BACKEND STREQUAL "OMP") 
    list(APPEND CUDA_NVCC_FLAGS
      -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP
      -Xcompiler -fopenmp)
    list(APPEND cuda_deps gomp)
  elseif (THRUST_BACKEND STREQUAL "TBB") 
    
  elseif (THRUST_BACKEND STREQUAL "CUDA")
    list(APPEND CUDA_NVCC_FLAGS
      -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA)
  endif ()

    #Eventually add some code instrumentation
  if (USE_NVTX)
    add_definitions(-DUSE_NVTX)
  endif ()

  #generic compilation flag, defines target architecture, type of build, ...
  list(APPEND CUDA_NVCC_FLAGS
    -gencode arch=compute_30,code=sm_30
    -gencode arch=compute_35,code=sm_35
    -gencode arch=compute_50,code=sm_50
    -gencode arch=compute_52,code=sm_52
    -gencode arch=compute_60,code=sm_60
    -gencode arch=compute_61,code=sm_61
    -std=c++11
    --expt-extended-lambda
    -rdc=true
    -Xcompiler -fPIC)

  set(CUDA_SEPARABLE_COMPILATION ON)

  if (CMAKE_BUILD_TYPE STREQUAL "Debug") 
    list(APPEND CUDA_NVCC_FLAGS 
      --device-debug
      --generate-line-info
      -g
      -G
      -DTHRUST_DEBUG )
  else ()
    list(APPEND CUDA_NVCC_FLAGS
      -use_fast_math
      -O3)
  endif ()
  #runtime libs for cuda
  list(APPEND cuda_deps cudadevrt nvToolsExt)
  link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64)
endif ()

############
# MPI Part #
############
find_package(MPI)
if (MPI_FOUND AND USE_MPI)
  add_definitions(-DUSE_MPI)
  include_directories(${MPI_INCLUDE_PATH})
  #runtime libs for mpi
  set(mpi_deps ${MPI_CXX_LIBRARIES})
endif ()
#############
# /MPI Part #
#############

###############
# Boost part  #
###############
find_package(Boost COMPONENTS mpi serialization) 
if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  link_directories(${Boost_LIBRARY_DIRS})
  #runtime libs for boost
  set(boost_deps boost_mpi boost_serialization)
endif ()
################
# /Boost part  #
################

#############
#  Testing  #
#############
if (TESTING_ENABLED)
  enable_testing()
endif ()

########################################
####        Doxygen Part            ####
########################################

# Add a target to generate API documentation with Doxygen
find_package(Doxygen)
if (DOXYGEN_FOUND AND USE_DOC)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(
    doc
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DOC_INSTALL_DIR}
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif ()

##############
#  Cppcheck  #
##############

# Add a target to generate Cppcheck report
include(${PRJ_CMAKE_DIR}/gen_cppcheck.cmake)
GENERATE_CPPCHECK(SOURCEDIRS src app test 
  ENABLE_IDS warning style performance portability information
  INCLUDES ${CURRENT_SOURCE_DIR}
  PROJECT_NAME ${PROJECT_NAME}
)

#####################################
#  Adding various code directories  #
#####################################

add_subdirectory(app)
add_subdirectory(binding)
add_subdirectory(src)

if (TESTING_ENABLED)
  add_subdirectory(test EXCLUDE_FROM_ALL)
endif ()
