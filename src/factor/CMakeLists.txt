add_library(factor
    primes_table.c
    trial_division.c)

target_include_directories(factor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
