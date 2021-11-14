cmake_minimum_required (VERSION 3.8)

find_package(Eigen3 CONFIG REQUIRED)

add_executable (bgv "bgv.cpp" "random_tools.cpp" "zq_type.h" "random_tools.h" "encryption.h")
target_link_libraries(bgv Eigen3::Eigen)
