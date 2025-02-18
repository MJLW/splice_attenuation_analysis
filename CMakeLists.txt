cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_STANDARD 99)

project(SpliceAttenuationAnalysis C)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/externals/sanitizers-cmake/cmake" ${CMAKE_MODULE_PATH})

include_directories(include)

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)
pkg_check_modules(CPLICEAI REQUIRED cpliceai)

add_subdirectory(src)

