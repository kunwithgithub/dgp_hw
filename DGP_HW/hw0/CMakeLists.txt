cmake_minimum_required(VERSION 3.1)
project(hw0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL        "Use OpenGL"                       ON)
option(LIBIGL_WITH_OPENGL_GLFW   "Use GLFW"                         ON)
option(LIBIGL_USE_STATIC_LIBRARY "Build LibIGL as a static library" OFF)

find_package(LIBIGL REQUIRED QUIET)

# Add your project files
add_executable(${PROJECT_NAME} main.cc Viewer.hh)
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw)
