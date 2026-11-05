
cmake_minimum_required(VERSION 3.20)

project(ENGO419_Lab8 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# Eigen 5.0.0
FetchContent_Declare(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/5.0.0/eigen-5.0.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

FetchContent_MakeAvailable(eigen)

# Boost 1.92.0
set(BOOST_INCLUDE_LIBRARIES math geometry)

FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.92.0/boost-1.92.0-cmake.tar.xz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

FetchContent_MakeAvailable(Boost)

# Executable
add_executable(ENGO419_Lab8
    main.cpp
    functions.cpp
    functions.h
)

target_link_libraries(ENGO419_Lab8 PRIVATE
    Eigen3::Eigen
    Boost::math
    Boost::geometry
)
