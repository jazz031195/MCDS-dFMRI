cmake_minimum_required(VERSION 3.1.0)
project(MC-DC_Simulator VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Git REQUIRED)
include(ExternalProject)

ExternalProject_Add(
    doctest
    URL https://github.com/doctest/doctest/archive/refs/tags/v2.4.11.tar.gz
    PREFIX ${CMAKE_BINARY_DIR}/externals/doctest
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/doctest"
)
set(DOCTEST_INCLUDE_DIR ${CMAKE_BINARY_DIR}/doctest/include/doctest)

ExternalProject_Add(
    eigen
    # 3.4.0 is not compatible with intellisense
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz
    PREFIX ${CMAKE_BINARY_DIR}/externals/eigen
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    PATCH_COMMAND ${GIT_EXECUTABLE} apply --ignore-whitespace -v ${CMAKE_CURRENT_SOURCE_DIR}/eigen_warnings.patch
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/eigen"
)
set(EIGEN_INCLUDE_DIR ${CMAKE_BINARY_DIR}/eigen/include/eigen3)


set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g0")
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3 -fsanitize=address -fsanitize=undefined")

set(SRC_DIR       MCDC_Simulator_public-master/src)
set(MAIN_TARGET   MC-DC_Simulator)
if (MSVC)  # Microsoft's cpp compiler
# TODO: To Be Tested
    set(COMPILER_OPTIONS
        ${COMPILER_ARCH_OPTIONS}
        -W4 -WX
    )
else()
    set(COMPILER_OPTIONS
        ${COMPILER_ARCH_OPTIONS}
        -Wall
        -Wextra
        # TODO: To make it fail on warning
        # -Werror
        -pedantic
        -Wwrite-strings
        -Winit-self
        -Wcast-qual
        -Wpointer-arith
        -Wstrict-aliasing
        -Wformat=2
        -Wuninitialized
        # Remove some of them
        -Wno-unused-parameter
    )
endif()
file(GLOB CPP_FILES ${SRC_DIR}/*.cpp)
list(REMOVE_ITEM CPP_FILES ${CMAKE_SOURCE_DIR}/${SRC_DIR}/main.cpp)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
add_executable(${MAIN_TARGET} ${CPP_FILES} ${SRC_DIR}/main.cpp)
# Compile release version with same flags as the test version (don't remove assertions)
# and only swap main() entry point. This is to avoid compiler warnings about unused variables in asserts
# Possible woraround:
# https://gitlab.mel.vin/mirror/doctest/-/commit/6b61e8aa3818c5ea100cedc1bb48a60ea10df6e8
target_compile_options(${MAIN_TARGET} PUBLIC ${COMPILER_OPTIONS} -DDONT_RUN_TESTS)
target_include_directories(${MAIN_TARGET} PUBLIC ${SRC_DIR} ${DOCTEST_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries(${MAIN_TARGET} PRIVATE Threads::Threads)
add_dependencies(${MAIN_TARGET} eigen doctest)

# Using a central precompiled header can increase compilation speed
# target_precompile_headers(${TARGET_APP_WITH_ETL} PUBLIC Source/main_header.hpp)
# TODO:
# if (TESTS)
#     # enable_testing()
#     # add_test(NAME ${MAIN_TARGET} COMMAND ${MAIN_TARGET})
#     # Make test executable
#     add_executable(tests tests.cpp)
#     target_compile_features(tests PRIVATE cxx_std_17)
#     target_link_libraries(tests PRIVATE doctest::doctest)
# endif()

add_executable(tests ${CPP_FILES})
target_compile_options(tests PUBLIC ${COMPILER_OPTIONS})
target_include_directories(tests PUBLIC ${SRC_DIR} ${DOCTEST_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries(tests Threads::Threads)
add_dependencies(tests eigen doctest)
