add_library(factor
    primes_table.c
    trial_division.c)

target_include_directories(factor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(${factor_include_dir} ${CMAKE_CURRENT_SOURCE_DIR})
