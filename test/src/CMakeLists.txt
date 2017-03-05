cmake_minimum_required(VERSION 3.5)

####################
#  Lib Dependancy  #
####################

include_directories( "${PRJ_ROOT_DIR}/src" )
set(lib_deps pdwt)

################
#  Build test  #
################
#add_lib_executable(libtest1 ${lib_deps})

if (CUDA_FOUND AND USE_CUDA)
  #add_cuda_lib_executable(cudalibtest1 ${lib_deps})
endif ()

##############
#  Unit test #
##############
#add_test(NAME libtest1 COMMAND libtest1)

if (CUDA_FOUND AND USE_CUDA)
#  add_test(NAME cudalibtest1 COMMAND cudalibtest1)
endif ()

##################
#  Config files  #
##################
set(DATA_FILES
  "${TEST_ROOT_DIR}/common/data/lena.dat"
)

install(FILES ${DATA_FILES}
  DESTINATION ${TEST_BIN_DIR} COMPONENT data)
