cmake_minimum_required (VERSION 3.5)

#####################################
#  Adding various code directories  #
#####################################

add_subdirectory(democpp)
add_subdirectory(demo2Dcpp)
add_subdirectory(demo3Dcpp)

if (CUDA_FOUND AND USE_CUDA)
  add_subdirectory(democuda)
endif ()
