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
prepare_executable(elliptic_curve elliptic.c)
prepare_executable(factor_trial_division factor_trial_division.c)
prepare_executable(fermat fermat.c)
prepare_executable(pollard_rho pollard_rho.c)
prepare_executable(pollard_rho2 pollard_rho2.c)
prepare_executable(pollard_rho3 pollard_rho3.c)
prepare_executable(pollard_p_1 pollard_p_1.c)
prepare_executable(shanks shanks.c)

#
# NUMBERS FOR THE TESTS
#

set(tiny1 1147 31 37)
set(tiny2 1679 23 73)
set(tiny3 137703491 17389 7919)

# both factors are bigger than the threshold set for the prime-table
set(tiny4 32588192363 180511 180533)

# The factors are the 33333th and the 40000th primes
set(tiny5 188695899619 393191 479909)

# 128 bits = 64 bits * 64 bits
set(128bit1 1004683540227020433292945938504490744493
    31696752165516604139 31696734573330559687)

# 108 bits = 64 bits * 44 bits
set(108bit1 408077400961771295988849803666327
    23196514630490635249 17592186044423)
set(108bit2 398897411515552667992550960001961
    21357450167740940629 18677202024709)
set(108bit3 1078836317125700231038586227595339
    33183360053317144537 32511364593347)
set(108bit4 576021769329963795775916997929453
    20043011670448667359 28739282239667)

# The number that shanks was trying to factor when he developed his
# SQUFOF algorithm, 2^60 + 2^30 - 1, the mean of M_31 and M_61
# the Mersenne primes
set(shanks1 1152921505680588799
    139001459 8294312261)

# Product of two twin primes: 2^50 + 2^30 + 2^10 ± 1
set(big1 1267653018083327397467897462783
    1125900980585471 1125900980585473)

do_test(prime_10 prime_test 10 29)
do_test(prime_15 prime_test 15 47)

do_test_factor(trial_division_tiny1 factor_trial_division ${tiny1})
do_test_factor(trial_division_tiny2 factor_trial_division ${tiny2})
do_test_factor(trial_division_tiny3 factor_trial_division ${tiny3})
do_test_factor(trial_division_tiny4 factor_trial_division ${tiny4})

do_test_factor(fermat_tiny4   fermat ${tiny4})
do_test_factor(fermat_tiny5   fermat ${tiny5})
do_test_factor(fermat_128bit1 fermat ${128bit1})

do_test_factor(pollard_rho_tiny5   pollard_rho ${tiny5})
do_test_factor(pollard_rho_108bit1 pollard_rho ${108bit1})

do_test_factor(pollard_rho2_tiny5   pollard_rho2 ${tiny5})
do_test_factor(pollard_rho2_108bit1 pollard_rho2 ${108bit1})

do_test_factor(pollard_rho3_tiny5   pollard_rho3 ${tiny5})
do_test_factor(pollard_rho3_108bit1 pollard_rho3 ${108bit1})

do_test_factor(pollard_p_1_tiny5   pollard_p_1 ${tiny5})
do_test_factor(pollard_p_1_108bit2 pollard_p_1 ${108bit2})
do_test_factor(pollard_p_1_108bit3 pollard_p_1 ${108bit3})
do_test_factor(pollard_p_1_108bit4 pollard_p_1 ${108bit4})

# Hard to find some number that the crappy algorithm for squfof I have works...
do_test_factor(shanks_tiny4  shanks  ${tiny4})

# Sill very slow, so don't do many tests...
do_test_factor(elliptic_curve_tiny5   elliptic_curve ${tiny5})
do_test_factor(elliptic_curve_shanks1 elliptic_curve ${shanks1})
