add_executable(prime_test prime_test.c)
add_executable(factor_trial_division factor_trial_division.c)
add_executable(fermat fermat.c)

target_link_libraries(prime_test factor m)
target_link_libraries(factor_trial_division factor m gmp)
target_link_libraries(fermat factor m gmp)
target_include_directories(prime_test PUBLIC ${factor_include_dir})
target_include_directories(factor_trial_division PUBLIC ${factor_include_dir})
target_include_directories(fermat PUBLIC ${factor_include_dir})

macro (do_test test_name test_executable arg result)
    add_test (${test_name} ${test_executable} ${arg})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION ${result})
endmacro (do_test)

do_test(prime_10 prime_test 10 29)
do_test(prime_15 prime_test 15 47)

# Two little tests to 'warm up the machine'
do_test(trial_division_1147 factor_trial_division 1147 31|37)
do_test(trial_division_1679 factor_trial_division 1679 23|73)

# Two bigger tests to see if everything works
# The second one uses prime numbers bigger than the higher threshold
# we have set for the default in the prime-table
do_test(trial_division_137703491 factor_trial_division 137703491 17389|7919)
do_test(trial_division_32588192363 factor_trial_division 32588192363 180511|180533)

do_test(fermat_32588192363 fermat 32588192363 180511|180533)
do_test(fermat_188695899619 fermat 188695899619 393191|479909)
