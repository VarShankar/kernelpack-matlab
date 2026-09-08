cmake_minimum_required(VERSION 3.18)

project(kernelpack_ibamr_rbc LANGUAGES C CXX)

find_package(IBAMR 0.19 REQUIRED)

add_executable(kernelpack_rbc_3d
    main.cpp
    RBCSurfaceMethod.cpp)

target_compile_features(kernelpack_rbc_3d PRIVATE cxx_std_17)
target_link_libraries(kernelpack_rbc_3d PRIVATE IBAMR::IBAMR3d)
