add_library(factor
    fermat.c
    pollard_rho.c
    pollard_p_1.c
    primes_table.c
    shanks.c
    trial_division.c)

target_include_directories(factor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(${factor_include_dir} ${CMAKE_CURRENT_SOURCE_DIR})
