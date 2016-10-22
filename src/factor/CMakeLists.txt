add_library(factor
    elliptic.c
    fermat.c
    pollard_rho.c
    pollard_p_1.c
    quadratic_sieve.c
    shanks.c
    trial_division.c)

target_include_directories(factor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(factor PUBLIC ${elliptic_include_dir})
target_include_directories(factor PUBLIC ${prime_include_dir})
target_link_libraries(factor elliptic prime)

set(${factor_include_dir} ${CMAKE_CURRENT_SOURCE_DIR})
