cmake_minimum_required(VERSION 3.1)
project(hw1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL" ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"   ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"  ON)
option(LIBIGL_USE_STATIC_LIBRARY "Build LibIGL as a static library" OFF)

find_package(LIBIGL REQUIRED QUIET)

# Build mode that compiles with optimizations but leaves in assertions.
set(CMAKE_CXX_FLAGS_RELWITHASSERT        "-O3" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
set(CMAKE_C_FLAGS_RELWITHASSERT          "-O3" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_RELWITHASSERT    "" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_RELWITHASSERT "" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
set(CMAKE_MODULE_LINKER_FLAGS_RELWITHASSERT "" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
set(CMAKE_STATIC_LINKER_FLAGS_RELWITHASSERT "" CACHE STRING "Flags used during RelWithAssertions builds" FORCE)
mark_as_advanced(
    CMAKE_CXX_FLAGS_RELWITHASSERT
    CMAKE_C_FLAGS_RELWITHASSERT
    CMAKE_EXE_LINKER_FLAGS_RELWITHASSERT
    CMAKE_SHARED_LINKER_FLAGS_RELWITHASSERT
    CMAKE_MODULE_LINKER_FLAGS_RELWITHASSERT
    CMAKE_STATIC_LINKER_FLAGS_RELWITHASSERT
)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to RelWithAssert")
    set(CMAKE_BUILD_TYPE "RelWithAssert")
endif()

# Make sure warnings/errors are still colorized when using Ninja for building.
add_definitions(-fdiagnostics-color=always)

# Add your project files
if (NOT SOLUTION)
    add_executable(${PROJECT_NAME} main.cc Viewer.hh extract_contours.cc extract_contours.hh)
    target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::opengl_glfw_imgui)
else()
    add_executable(${PROJECT_NAME} main_solution.cc Viewer.hh extract_contours_solution.cc extract_contours.hh)
    target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::opengl_glfw_imgui)
endif()
