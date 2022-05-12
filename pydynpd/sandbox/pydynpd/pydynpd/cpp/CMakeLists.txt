cmake_minimum_required(VERSION 3.1)
##cd build
##cmake ..
##make install

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE)
project(gmm_module VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "-Wall -O3 -m64 ")
#-m64 -fopenmp -march=native  -mavx2 
execute_process (
    COMMAND bash -c "python3 -m pybind11 --includes"
    OUTPUT_VARIABLE outVar
)

# Include dir
include_directories(/usr/local/include)
include_directories(/usr/include/eigen3)
include_directories(${outVar})
include_directories(/usr/include/python3.9)
include_directories(/usr/local/lib/python3.9/dist-packages/pybind11/include)
include_directories(/usr/include/python3.10)
include_directories(/usr/local/lib/python3.10/dist-packages/pybind11/include)




#find_package(pybind11 REQUIRED)



# Headers
#set(PROJECT_SOURCE_DIR "/src")
#set(PROJECT_INCLUDE_DIR "/src")

set(LIB_OUT "lib")

# Source files
set(SOURCE_FILES

    Command.cpp GMM.cpp Specification_Test.cpp  dynamic_model.cpp   info.cpp  variable.cpp Common_Functions.cpp   List_Variables.cpp  Step_Result.cpp         instruments.cpp    wrap.cpp
)

add_library(gmm_module SHARED  ${SOURCE_FILES})


# Command.cpp GMM.cpp Specification_Test.cpp  dynamic_model.cpp  variable.cpp Common_Functions.cpp   
# List_Variables.cpp  Step_Result.cpp         instruments.cpp    wrap.cpp

 set_target_properties(gmm_module PROPERTIES PREFIX "")

set_target_properties(gmm_module PROPERTIES OUTPUT_NAME "gmm_module")

install(TARGETS gmm_module DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../)