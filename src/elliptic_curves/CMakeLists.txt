add_library(elliptic
    elliptic.c)

target_include_directories(elliptic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(${elliptic_include_dir} ${CMAKE_CURRENT_SOURCE_DIR})
