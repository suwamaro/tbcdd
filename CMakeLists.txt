cmake_minimum_required(VERSION 2.8)

project(tbc)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


# -- COMPILE OPTIONS -------------------------------------------

# Performance/debug options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -gdwarf-3 -Wall -DADD_")
add_definitions(-D_GLIBCXX_DEBUG=1) # libstdc++
add_definitions(-D_LIBCPP_DEBUG2=1) # libc++

# C++17
include(CheckCXXCompilerFlag)
#check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORTS_CXX17)
#if(NOT COMPILER_SUPPORTS_CXX17)
#  message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++17 support.# Aborting.")
#endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

# Mac fixes
if (APPLE)
  cmake_policy(SET CMP0042 NEW) # @rpath
  set(LIBS ${LIBS} "-framework Accelerate")
  set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
endif()


# -- DEPENDENCIES -------------------------------------------

# Armadillo
find_package(Armadillo)
include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})
set(LIBDIR ${LIBDIR} ${ARMADILLO_LIBRARY_DIRS})

# CPPTOML
find_package(CPPTOML REQUIRED)
include_directories(SYSTEM ${CPPTOML_INCLUDE_DIRS})


# -- BUILD -------------------------------------------

include_directories(src)
link_directories(${LIBDIR})
add_subdirectory(src)

# -- REPORT -------------------------------------------

message("")
message("Configuration results")
message("---------------------")
message("C++ compiler        : ${CMAKE_CXX_COMPILER}")
message("C++ compiler flags  : ${CMAKE_CXX_FLAGS}")

get_property(INC_DIRS DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("LIBS:     ${LIBS}")
message("INCLUDES: ${INC_DIRS}")
message("WITH_MPI: ${MPI_FOUND}")
message("")
