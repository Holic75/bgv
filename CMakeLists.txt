cmake_minimum_required (VERSION 3.8)
project ("BGV")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

# Include sub-projects.
add_subdirectory("src")

# TODO: Add tests and install targets if needed.
