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
