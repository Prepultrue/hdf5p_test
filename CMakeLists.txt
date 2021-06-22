cmake_minimum_required(VERSION 3.19)
project(testc C)

set(CMAKE_C_STANDARD 11)

find_package(HDF5 REQUIRED)
find_package(MPI REQUIRED)

include_directories(${HDF5_INCLUDE_DIRS} ${MPI_INCLUDE_DIRS})

add_executable(testc main.c)
target_link_libraries(testc ${HDF5_LIBRARIES} ${MPI_LIBRARIES})
