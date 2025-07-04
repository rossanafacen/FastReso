CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(fastreso CXX)

set(CMAKE_MODULE_PATH ${fastreso_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -fopenmp -Wall -std=c++20")

#------------------------------------------------------------------------
# Turn on GSL
#------------------------------------------------------------------------
find_package(GSL)
include_directories(${GSL_INCLUDE_DIRS})

##------------------------------------------------------------------------
## Add source folder
##------------------------------------------------------------------------

include_directories(./)


# defines target library dist
set(fastreso_source
    TParticle.cpp
    TFastReso_THERMINATOR.cpp
    TFastReso_AZYHYDRO.cpp
    TParticle_THERMINATOR.cpp
    TParticle_AZYHYDRO.cpp
    TKernel.cpp
    )

add_library(fastreso ${fastreso_source})
include_directories("/cvmfs/vae.gsi.de/vae25/spack-0.23/opt/linux-debian12-zen/gcc-13.3.0/nlohmann-json-3.11.3-wbq27t5n65quivp5xfbkehtxcdi4xfrp/include")

target_link_libraries(fastreso ${GSL_LIBRARIES})

option(FASTRESO_ROOT "Link with ROOT Librariers" OFF)
if (FASTRESO_ROOT)
    message("FASTRESO_ROOT is on") 
   add_definitions("-DEKT_ROOT") 
   find_package(ROOT)
   include_directories(AFTER ${ROOT_INCLUDE_DIR})
   message("FASTRESO_ROOT links with ${ROOT_LIBRARIES}")
endif()



# compile executables
add_executable(run_fastreso_therminator.exe run_fastreso_therminator.cpp)
add_executable(run_fastreso_azyhydro.exe run_fastreso_azyhydro.cpp)
add_executable(run_kernels.exe run_kernels.cpp)
add_executable(run_kernels_PCE.exe run_kernels_PCE.cpp)
add_executable(run_eos.exe run_eos.cpp)
add_executable(run_PCE.exe run_PCE.cpp)

# link libraries
target_link_libraries(run_fastreso_therminator.exe fastreso)
target_link_libraries(run_fastreso_azyhydro.exe fastreso)
target_link_libraries(run_kernels.exe fastreso)
target_link_libraries(run_kernels_PCE.exe fastreso)
target_link_libraries(run_eos.exe fastreso)
target_link_libraries(run_PCE.exe fastreso)


if (FASTRESO_ROOT)
add_executable(run_kernels_root.exe run_kernels_root.cpp)
target_link_libraries(run_kernels_root.exe fastreso)
target_link_libraries(run_kernels_root.exe ${ROOT_LIBRARIES})
endif()




