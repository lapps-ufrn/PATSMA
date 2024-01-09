add_executable (integrity 
  main-test.cpp 
  test-integrity.cpp
)
target_include_directories (integrity PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_link_libraries (integrity autotuning)
install (TARGETS integrity DESTINATION bin)
