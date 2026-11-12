cmake_minimum_required(VERSION 3.20)
project(TextureConsolidation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extern/Misha/Misha/Geometry.h OR
   NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extern/PhongRodriguesVF/Include/SimplicialMeshProcessing.h OR
   NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extern/TextureSignalProcessing/TextureStitching/TextureStitching.cpp)
  message(FATAL_ERROR "Submodules missing. Run: git submodule update --init")
endif()

find_package(OpenMP REQUIRED)
find_package(PNG REQUIRED)
find_package(JPEG REQUIRED)
find_package(ZLIB REQUIRED)

# Compiler settings and libraries shared by every executable.
add_library(tc_settings INTERFACE)
target_link_libraries(tc_settings INTERFACE OpenMP::OpenMP_CXX PNG::PNG JPEG::JPEG ZLIB::ZLIB)
# The vendored Eigen 3.3.7 (in PhongRodriguesVF and TextureSignalProcessing) uses
# std::result_of when it believes the compiler has it, and C++20 removed it; libc++ (Apple
# clang) enforces that. Eigen's own fallback is selected by this macro, which it guards
# with #ifndef for exactly this purpose.
target_compile_definitions(tc_settings INTERFACE EIGEN_HAS_STD_RESULT_OF=0)
if(WIN32)
  # Misha includes <png.h>/<jpeglib.h>/<zlib.h> on Windows only when USE_VCPKG is set.
  target_compile_definitions(tc_settings INTERFACE USE_VCPKG _CRT_SECURE_NO_WARNINGS)
endif()
if(MSVC)
  target_compile_options(tc_settings INTERFACE /bigobj /permissive- /W1 /wd4819) # wd4819: non-ASCII characters in comments of a Misha header
else()
  target_compile_options(tc_settings INTERFACE -Wno-deprecated -Wno-unused-result)
endif()

# Our programs. Include order matters: extern/Misha must precede PhongRodriguesVF's
# vendored copy of Misha (third_party/include), which is then only used for Eigen.
add_library(tc_common INTERFACE)
target_include_directories(tc_common INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/Misha
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/PhongRodriguesVF
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/PhongRodriguesVF/third_party/include)
target_link_libraries(tc_common INTERFACE tc_settings)

function(tc_add_executable name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} PRIVATE tc_common)
endfunction()

tc_add_executable(BackProjection src/BackProjection/BackProjection.cpp)
tc_add_executable(TextureFlow src/TextureFlow/TextureFlow.cpp)

# Gradient-domain stitching (Sec. 6.3, Eq. 7), built unmodified from TextureSignalProcessing (without its OpenGL viewer).
# extern/Misha precedes TSP's vendored copy of Misha, which is then only used for Eigen.
set(TSP ${CMAKE_CURRENT_SOURCE_DIR}/extern/TextureSignalProcessing)
# Built through a one-line wrapper that makes Misha's MK_ERROR_OUT throw, so that every
# failure ends with a non-zero exit status (see the wrapper for why).
add_executable(TextureDomainStitching src/TextureDomainStitching/TextureDomainStitching.cpp)
target_include_directories(TextureDomainStitching PRIVATE
  ${TSP}
  ${TSP}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/Misha
  ${TSP}/third_party/include)
target_compile_definitions(TextureDomainStitching PRIVATE NO_OPEN_GL_VISUALIZATION)
target_link_libraries(TextureDomainStitching PRIVATE tc_settings)

