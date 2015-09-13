macro(prepare_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} montgomery_multiplication m gmp)
    target_include_directories(${target} PUBLIC ${montgomery_multiplication_include_dir})
endmacro(prepare_executable)

macro(do_test_mul test_name test_executable
        op1 op2 mod result)
    add_test (${test_name} ${test_executable} ${op1} ${op2} ${mod})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION ${result})
endmacro(do_test_mul)

prepare_executable(mont_mul mont_mul.c)
prepare_executable(mont_pow mont_pow.c)
prepare_executable(mont_pow_mpz mont_pow_mpz.c)

#
# NUMBERS FOR THE TESTS
#

set(mont_mul_numbers1 12345 54321 55555 43895)
set(mont_mul_numbers2 11111 99999 55779 26988)

set(mont_mul_numbers3 123456789 987654321 1234567654321 556733621704)
set(mont_mul_numbers4 9876789876789 9977553311 1234567654321 770820326330)

do_test_mul(mont_mul1 mont_mul ${mont_mul_numbers1})
do_test_mul(mont_mul2 mont_mul ${mont_mul_numbers2})
do_test_mul(mont_mul3 mont_mul ${mont_mul_numbers3})
do_test_mul(mont_mul4 mont_mul ${mont_mul_numbers4})

set(mont_pow_numbers1_1 12345 2 55555 11660)
set(mont_pow_numbers1_2 12345 3 55555 55250)
set(mont_pow_numbers1_3 12345 4 55555 12515)
set(mont_pow_numbers1_4 12345 5 55555 54775)
set(mont_pow_numbers1_5 12345 6 55555 37470)

set(mont_pow_numbers2_1 9876789876 2 1234567654321 233990428508)
set(mont_pow_numbers2_2 9876789876 3 1234567654321 793520975495)
set(mont_pow_numbers2_3 9876789876 4 1234567654321 947178751242)
set(mont_pow_numbers2_4 9876789876 5 1234567654321 441750524501)
set(mont_pow_numbers2_5 9876789876 6 1234567654321 477805590232)

do_test_mul(mont_pow1_1 mont_pow ${mont_pow_numbers1_1})
do_test_mul(mont_pow1_2 mont_pow ${mont_pow_numbers1_2})
do_test_mul(mont_pow1_3 mont_pow ${mont_pow_numbers1_3})
do_test_mul(mont_pow1_4 mont_pow ${mont_pow_numbers1_4})
do_test_mul(mont_pow1_5 mont_pow ${mont_pow_numbers1_5})

do_test_mul(mont_pow2_1 mont_pow ${mont_pow_numbers2_1})
do_test_mul(mont_pow2_2 mont_pow ${mont_pow_numbers2_2})
do_test_mul(mont_pow2_3 mont_pow ${mont_pow_numbers2_3})
do_test_mul(mont_pow2_4 mont_pow ${mont_pow_numbers2_4})
do_test_mul(mont_pow2_5 mont_pow ${mont_pow_numbers2_5})

# exponents are 2^64 + 2^32 - 1 and 2^66 + 2^33 - 1
set(mont_pow_mpz_numbers2_1 9876789876 18446744078004518911 1234567654321 105687510110)
set(mont_pow_mpz_numbers2_2 9876789876 73786976303428141055 1234567654321 685776966138)

do_test_mul(mont_pow_mpz1_1 mont_pow_mpz ${mont_pow_numbers1_1})
do_test_mul(mont_pow_mpz1_2 mont_pow_mpz ${mont_pow_numbers1_2})
do_test_mul(mont_pow_mpz1_3 mont_pow_mpz ${mont_pow_numbers1_3})
do_test_mul(mont_pow_mpz1_4 mont_pow_mpz ${mont_pow_numbers1_4})
do_test_mul(mont_pow_mpz1_5 mont_pow_mpz ${mont_pow_numbers1_5})

do_test_mul(mont_pow_mpz2_1 mont_pow_mpz ${mont_pow_mpz_numbers2_1})
do_test_mul(mont_pow_mpz2_2 mont_pow_mpz ${mont_pow_mpz_numbers2_2})
