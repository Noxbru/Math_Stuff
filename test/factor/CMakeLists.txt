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

macro(do_test_factor test_name test_executable arg result1 result2)
    do_test(${test_name} ${test_executable} ${arg} ${result1}|${result2})
endmacro(do_test_factor)

prepare_executable(prime_test prime_test.c)
prepare_executable(factor_trial_division factor_trial_division.c)
prepare_executable(fermat fermat.c)
prepare_executable(pollard_rho pollard_rho.c)
prepare_executable(pollard_rho2 pollard_rho2.c)
prepare_executable(pollard_p_1 pollard_p_1.c)

do_test(prime_10 prime_test 10 29)
do_test(prime_15 prime_test 15 47)

# Two little tests to 'warm up the machine'
do_test_factor(trial_division_1147 factor_trial_division 1147 31 37)
do_test_factor(trial_division_1679 factor_trial_division 1679 23 73)

# Two bigger tests to see if everything works
# The second one uses prime numbers bigger than the higher threshold
# we have set for the default in the prime-table
do_test_factor(trial_division_137703491   factor_trial_division 137703491    17389   7919)
do_test_factor(trial_division_32588192363 factor_trial_division 32588192363 180511 180533)

do_test_factor(fermat_32588192363  fermat  32588192363 180511 180533)
do_test_factor(fermat_188695899619 fermat 188695899619 393191 479909)

do_test_factor(fermat_64bit fermat
    1004683540227020433292945938504490744493
    31696752165516604139
    31696734573330559687)

do_test_factor(pollard_rho_188695899619 pollard_rho 188695899619 393191 479909)
do_test_factor(pollard_rho_108bit pollard_rho
    408077400961771295988849803666327   # 108 bit number
    23196514630490635249                #  64 bit number
    17592186044423)                     #  44 bit number

do_test_factor(pollard_rho2_188695899619 pollard_rho2 188695899619 393191 479909)
do_test_factor(pollard_rho2_108bit pollard_rho2
    408077400961771295988849803666327   # 108 bit number
    23196514630490635249                #  64 bit number
    17592186044423)                     #  44 bit number

do_test_factor(pollard_p_1_188695899619 pollard_p_1 188695899619 393191 479909)
do_test_factor(pollard_p_1_108bit pollard_p_1
    398897411515552667992550960001961   # 108 bit number
    21357450167740940629                #  64 bit number
    18677202024709)                     #  44 bit number
