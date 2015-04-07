add_executable(prime_test prime_test.c)

target_link_libraries(prime_test factor m)
target_include_directories(prime_test PUBLIC ${factor_include_dir})

macro (do_test test_name test_executable arg result)
    add_test (${test_name} ${test_executable} ${arg})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION ${result})
endmacro (do_test)

do_test(prime_10 prime_test 10 29)
do_test(prime_15 prime_test 15 47)
