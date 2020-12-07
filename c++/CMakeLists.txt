project(hw02)

cmake_minimum_required(VERSION 3.2)

add_definitions(-std=c++11)

include_directories( ${CMAKE_SOURCE_DIR}/thirdparty )

FILE(GLOB SRC_FILES src/*.cpp)
add_executable(hw02 ${SRC_FILES})

if (WIN32)
    set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT "hw02")
endif (WIN32)
