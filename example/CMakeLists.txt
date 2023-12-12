add_executable (rb_solve 
  parallel_rb.cpp
)

target_include_directories (rb_solve PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_link_libraries (rb_solve autotuning)
install (TARGETS rb_solve DESTINATION bin)


