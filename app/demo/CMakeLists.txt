cmake_minimum_required(VERSION 3.5)

##################
#  Generic Part  #
##################
project(demo)
set(OUTPUT_NAME  demo)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

###############
#  Src files  #
###############

set(sources
  demo.cpp
  io.cpp
)

set(headers
)

####################
#  Lib Dependancy  #
####################

include_directories( "../../src" )

################
#  Executable  #
################

add_executable(${OUTPUT_NAME} ${sources})
target_link_libraries(${OUTPUT_NAME}
  pdwt
  ${mpi_deps}
  ${boost_deps})
set_property(TARGET ${OUTPUT_NAME}
  PROPERTY OUTPUT_NAME ${OUTPUT_NAME})
install(TARGETS ${OUTPUT_NAME}
  RUNTIME DESTINATION ${RUNTIME_DEST})
