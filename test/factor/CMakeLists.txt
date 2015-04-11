macro(prepare_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} factor m gmp)
    target_include_directories(${target} PUBLIC ${factor_include_dir})
endmacro(prepare_executable)

macro (do_test test_name test_executable arg result)
    add_test (${test_name} ${test_executable} ${arg})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION ${result})
endmacro (do_test)

prepare_executable(prime_test prime_test.c)
prepare_executable(factor_trial_division factor_trial_division.c)
prepare_executable(fermat fermat.c)

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
