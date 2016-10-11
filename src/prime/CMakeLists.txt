add_library(prime
    primes_table.c)

target_include_directories(factor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(${prime_include_dir} ${CMAKE_CURRENT_SOURCE_DIR})
