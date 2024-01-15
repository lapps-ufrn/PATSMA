# add_executable (inside RBGaussSeidel-insideAT.cpp)
# target_include_directories (inside PUBLIC "${PROJECT_SOURCE_DIR}/src")
# target_link_libraries (inside autotuning)
# install (TARGETS inside DESTINATION bin)

# add_executable (outside RBGaussSeidel-outsideAT.cpp)
# target_include_directories (outside PUBLIC "${PROJECT_SOURCE_DIR}/src")
# target_link_libraries (outside autotuning)
# install (TARGETS outside DESTINATION bin)

# add_executable (no RBGaussSeidel-noAT.cpp)
# install (TARGETS no DESTINATION bin)

cmake_minimum_required (VERSION 3.5)

project(RBGaussExemple)

find_package(OpenMP QUIET)
if(NOT OPENMP_FOUND)
  message(FATAL_ERROR "OpenMP must exist to compile FWI when OpenMP is enabled.")
else()
  message(STATUS "OpenMP ${OpenMP_CXX_VERSION} found & enabled.")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

# Find and link the library created in the PATSMA project
find_package(PATSMA_lib REQUIRED HINTS patsma)
# Check if the library was found
if(NOT PATSMA_FOUND)
  message(FATAL_ERROR "PATSMA library not found")
endif()

# Create the example executable and link the library
set(run_name insideat)
add_executable(${run_name} RBGaussSeidel-insideAT.cpp)
target_include_directories (${run_name} PUBLIC ${PATSMA_INCLUDE_DIRS})
target_link_libraries (${run_name} ${PATSMA_LIBRARIES})
install (TARGETS ${run_name} DESTINATION bin)

set(run_name outsideat)
add_executable(${run_name} RBGaussSeidel-outsideAT.cpp)
target_include_directories (${run_name} PUBLIC ${PATSMA_INCLUDE_DIRS})
target_link_libraries (${run_name} ${PATSMA_LIBRARIES})
install (TARGETS ${run_name} DESTINATION bin)

set(run_name noat)
add_executable(${run_name} RBGaussSeidel-noAT.cpp)
target_include_directories (${run_name} PUBLIC ${PATSMA_INCLUDE_DIRS})
target_link_libraries (${run_name} ${PATSMA_LIBRARIES})
install (TARGETS ${run_name} DESTINATION bin)