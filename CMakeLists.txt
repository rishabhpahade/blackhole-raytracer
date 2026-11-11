cmake_minimum_required(VERSION 3.14)
project(blackhole_raytracer LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler optimization and warning flags
if(MSVC)
    add_compile_options(/W4 /O2 /fp:precise)
else()
    # NOTE: do not add -ffast-math here. It implies -ffinite-math-only, which
    # lets the compiler fold std::isfinite() to true and silently disables the
    # RayStatus::NumericalFailure guard in GeodesicIntegrator.
    add_compile_options(-Wall -Wextra -Wpedantic -O3 -fno-math-errno)
endif()

# Optional OpenMP for parallel rendering
find_package(OpenMP)

# Header include paths
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Core library
add_library(blackhole_core STATIC
    src/Camera.cpp
    src/AccretionDisk.cpp
    src/GeodesicIntegrator.cpp
    src/Renderer.cpp
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(blackhole_core PUBLIC OpenMP::OpenMP_CXX)
endif()

# Main Ray Tracer Executable
add_executable(blackhole_raytracer src/main.cpp)
target_link_libraries(blackhole_raytracer PRIVATE blackhole_core)

# Scientific Validation Test Suite
enable_testing()
add_executable(validation_tests tests/validation_tests.cpp)
target_link_libraries(validation_tests PRIVATE blackhole_core)
add_test(NAME ScientificValidation COMMAND validation_tests)
