add_library(montgomery_multiplication
    mont.c)

target_include_directories(montgomery_multiplication PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(${montgomery_multiplication_include_dir} ${CMAKE_CURRENT_SOURCE_DIR})
