# add_executable (integrity
#   main-test.cpp
#   test-integrity.cpp
# )
# target_include_directories (integrity PUBLIC "${PROJECT_SOURCE_DIR}/src")
# target_link_libraries (integrity autotuning)
# install (TARGETS integrity DESTINATION bin)

cmake_minimum_required (VERSION 3.5)

project(test)

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
set(run_name integrity)
add_executable(${run_name} main-test.cpp test-integrity.cpp)
target_include_directories (${run_name} PUBLIC ${PATSMA_INCLUDE_DIRS})
target_link_libraries (${run_name} ${PATSMA_LIBRARIES})
install (TARGETS ${run_name} DESTINATION bin)