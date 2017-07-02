cmake_minimum_required (VERSION 3.5)

#####################################
#  Adding binding code directories  #
#####################################

include_directories("pybind11")
add_subdirectory(pybind11)
add_subdirectory(src)
