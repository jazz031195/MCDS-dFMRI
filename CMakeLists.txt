cmake_minimum_required(VERSION 3.1.0)
project(MC-DC_Simulator VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
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
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
add_executable(${MAIN_TARGET} ${CPP_FILES})
target_compile_options(${MAIN_TARGET} PUBLIC ${COMPILER_OPTIONS})
target_include_directories(${MAIN_TARGET} PUBLIC ${SRC_DIR})
target_link_libraries(${MAIN_TARGET} PRIVATE Threads::Threads)
# Using a central precompiled header can increase compilation speed
# target_precompile_headers(${TARGET_APP_WITH_ETL} PUBLIC Source/main_header.hpp)
# TODO:
if (TESTS)
    enable_testing()
    add_test(NAME ${MAIN_TARGET} COMMAND ${MAIN_TARGET})
endif()