cmake_minimum_required(VERSION 3.16)

project(tests)

add_executable(${PROJECT_NAME} finite_differences.cpp)

# Link to DCA and make the header files accessible only to this executable.
target_link_libraries(${PROJECT_NAME} PRIVATE DCA::DCA)

# The following lines are only needed if an app needs to know if CompactNSearch is built.
# Normally, they can be neglected.
if(${USE_COMPACT_N_SEARCH})
  target_compile_definitions(${PROJECT_NAME} PRIVATE -DBUILD_COMPACT_N_SEARCH=1)
else(${USE_COMPACT_N_SEARCH})
  target_compile_definitions(${PROJECT_NAME} PRIVATE -DBUILD_COMPACT_N_SEARCH=0)
endif(${USE_COMPACT_N_SEARCH})