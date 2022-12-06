function(example example_file)
  add_executable(${example_file} ${example_file}.c)
  target_link_libraries(${example_file} wildriver)
endfunction()

example(size)
