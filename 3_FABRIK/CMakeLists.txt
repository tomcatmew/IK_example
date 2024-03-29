cmake_minimum_required(VERSION 3.12)

######################################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
IF (MSVC)
  set(CMAKE_CXX_FLAGS "/W4")
ELSE ()
  set(CMAKE_CXX_FLAGS "-Wall -g")
ENDIF ()

####################################

project(3_FABRIK)

add_definitions(-DPATH_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

# glfw
SET(GLFW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../3rd_party/libglfw)
include(${CMAKE_CURRENT_SOURCE_DIR}/../3rd_party/FindGLFW.cmake)

# opengl
find_package(OpenGL REQUIRED)

# delfem2
add_definitions(-DDFM2_HEADER_ONLY=ON)
set(DELFEM2_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../3rd_party/delfem2/include")

# imgui_glfw
add_definitions(-DIMGUI_IMPL_OPENGL_LOADER_GLAD)
set(IMGUI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../3rd_party/imgui")

####################################################

include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLFW_INCLUDE_DIR}
    ${DELFEM2_INCLUDE_DIR}
    ${IMGUI_DIR}
    ${IMGUI_DIR}/backends
)

add_executable(${PROJECT_NAME}
    ${IMGUI_DIR}/imgui.h
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl2.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl2.h
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.h
    main.cpp
    )

target_link_libraries(${PROJECT_NAME}
    ${DFM2_LIBRARY}
    ${OPENGL_LIBRARY}
    ${GLFW_LIBRARIES}
    )
