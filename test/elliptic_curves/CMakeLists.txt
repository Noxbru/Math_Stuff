macro(prepare_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} elliptic m gmp)
    target_include_directories(${target} PUBLIC ${elliptic_include_dir})
endmacro(prepare_executable)

macro(do_test_curve_sum test_name test_executable
        c_a c_b c_m p1_x p1_y p2_x p2_y r_x r_y)
    add_test (${test_name} ${test_executable}
        ${c_a} ${c_b} ${c_m} ${p1_x} ${p1_y} ${p2_x} ${p2_y})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION "${r_x}\t${r_y}")
endmacro(do_test_curve_sum)

macro(do_test_curve_double test_name test_executable
        c_a c_b c_m p1_x p1_y r_x r_y)
    add_test (${test_name} ${test_executable}
        ${c_a} ${c_b} ${c_m} ${p1_x} ${p1_y})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION "${r_x}\t${r_y}")
endmacro(do_test_curve_double)

macro(do_test_curve_mul test_name test_executable
        c_a c_b c_m p1_x p1_y times r_x r_y)
    add_test (${test_name} ${test_executable}
        ${c_a} ${c_b} ${c_m} ${p1_x} ${p1_y} ${times})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION "${r_x}\t${r_y}")
endmacro(do_test_curve_mul)

prepare_executable(elliptic_sum elliptic_sum.c)
prepare_executable(elliptic_double elliptic_double.c)
prepare_executable(elliptic_mul elliptic_mul.c)

set(curve1 5 -5 65537)
set(curve1_point1     1     1)
set(curve1_point2    14 65484) # 2 * (1,1)
set(curve1_point3 24821 62767) # 3 * (1,1)
set(curve1_point4 42740 22681) # 4 * (1,1)
set(curve1_point5 50687 11143) # 5 * (1,1)
set(curve1_point6 45879 59979) # 6 * (1,1)

do_test_curve_sum(curve_sum1 elliptic_sum ${curve1}
    ${curve1_point1} ${curve1_point1} ${curve1_point2})
do_test_curve_sum(curve_sum2 elliptic_sum ${curve1}
    ${curve1_point1} ${curve1_point2} ${curve1_point3})
do_test_curve_sum(curve_sum3 elliptic_sum ${curve1}
    ${curve1_point1} ${curve1_point3} ${curve1_point4})
do_test_curve_sum(curve_sum4 elliptic_sum ${curve1}
    ${curve1_point2} ${curve1_point3} ${curve1_point5})

do_test_curve_double(curve_double1 elliptic_double ${curve1}
    ${curve1_point1} ${curve1_point2})
do_test_curve_double(curve_double2 elliptic_double ${curve1}
    ${curve1_point2} ${curve1_point4})

do_test_curve_mul(curve_mul1 elliptic_mul ${curve1}
    ${curve1_point1} 3 ${curve1_point3})
do_test_curve_mul(curve_mul2 elliptic_mul ${curve1}
    ${curve1_point1} 5 ${curve1_point5})
do_test_curve_mul(curve_mul3 elliptic_mul ${curve1}
    ${curve1_point2} 3 ${curve1_point6})
