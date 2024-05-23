cmake_minimum_required(VERSION 3.28)
project(Coursework)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/data)

# Add the executable
add_executable(Coursework
        data/main.cpp
        data/csv.cpp
        data/read_data.cpp
        data/linalg.cpp
        data/Markowitz.cpp
        data/write_data.cpp
)

# Copy asset_returns.csv to the build directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/asset_returns.csv ${CMAKE_CURRENT_BINARY_DIR}/asset_returns.csv COPYONLY)
