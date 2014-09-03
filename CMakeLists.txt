cmake_minimum_required(VERSION 2.6)
project(libigl_tutorials)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

## Check for CoMiSo, if not available skip the examples that depends on it
find_package(LIBCOMISO QUIET)

## Check for MATLAB, if not available skip the examples that depends on it
find_package(MATLAB QUIET)

## Check for EMBREE, if not available skip the examples that depends on it
find_package(EMBREE QUIET)

add_subdirectory("demo")
