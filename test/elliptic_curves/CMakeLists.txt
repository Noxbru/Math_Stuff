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

prepare_executable(elliptic_sum_weierstrass_affine    elliptic_sum_weierstrass_affine.c)
prepare_executable(elliptic_double_weierstrass_affine elliptic_double_weierstrass_affine.c)
prepare_executable(elliptic_mul_weierstrass_affine    elliptic_mul_weierstrass_affine.c)

prepare_executable(elliptic_sum_montgomery_affine    elliptic_sum_montgomery_affine.c)
prepare_executable(elliptic_double_montgomery_affine elliptic_double_montgomery_affine.c)
prepare_executable(elliptic_mul_montgomery_affine    elliptic_mul_montgomery_affine.c)

set(weir_curve1 5 -5 65537)
set(weir_curve1_point1     1     1)
set(weir_curve1_point2    14 65484) # 2 * (1,1)
set(weir_curve1_point3 24821 62767) # 3 * (1,1)
set(weir_curve1_point4 42740 22681) # 4 * (1,1)
set(weir_curve1_point5 50687 11143) # 5 * (1,1)
set(weir_curve1_point6 45879 59979) # 6 * (1,1)

do_test_curve_sum(elliptic_curve_weierstrass_sum1 elliptic_sum_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} ${weir_curve1_point1} ${weir_curve1_point2})
do_test_curve_sum(elliptic_curve_weierstrass_sum2 elliptic_sum_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} ${weir_curve1_point2} ${weir_curve1_point3})
do_test_curve_sum(elliptic_curve_weierstrass_sum3 elliptic_sum_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} ${weir_curve1_point3} ${weir_curve1_point4})
do_test_curve_sum(elliptic_curve_weierstrass_sum4 elliptic_sum_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point2} ${weir_curve1_point3} ${weir_curve1_point5})

do_test_curve_double(elliptic_curve_weierstrass_double1 elliptic_double_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} ${weir_curve1_point2})
do_test_curve_double(elliptic_curve_weierstrass_double2 elliptic_double_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point2} ${weir_curve1_point4})

do_test_curve_mul(elliptic_curve_weierstrass_mul1 elliptic_mul_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} 3 ${weir_curve1_point3})
do_test_curve_mul(elliptic_curve_weierstrass_mul2 elliptic_mul_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} 5 ${weir_curve1_point5})
do_test_curve_mul(elliptic_curve_weierstrass_mul3 elliptic_mul_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point2} 3 ${weir_curve1_point6})

set(mont_curve1 5 3 65537)
set(mont_curve1_point1      3     5)
set(mont_curve1_point2   6117  1282) # 2 * (3,5)
set(mont_curve1_point3  54943 58216) # 3 * (3,5)
set(mont_curve1_point4  62979 59844) # 4 * (3,5)
set(mont_curve1_point5  23228 18695) # 5 * (3,5)
set(mont_curve1_point6  12508 22665) # 6 * (3,5)

do_test_curve_sum(elliptic_curve_montgomery_sum1 elliptic_sum_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} ${mont_curve1_point1} ${mont_curve1_point2})
do_test_curve_sum(elliptic_curve_montgomery_sum2 elliptic_sum_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} ${mont_curve1_point2} ${mont_curve1_point3})
do_test_curve_sum(elliptic_curve_montgomery_sum3 elliptic_sum_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} ${mont_curve1_point3} ${mont_curve1_point4})
do_test_curve_sum(elliptic_curve_montgomery_sum4 elliptic_sum_montgomery_affine ${mont_curve1}
    ${mont_curve1_point2} ${mont_curve1_point3} ${mont_curve1_point5})

do_test_curve_double(elliptic_curve_montgomery_double1 elliptic_double_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} ${mont_curve1_point2})
do_test_curve_double(elliptic_curve_montgomery_double2 elliptic_double_montgomery_affine ${mont_curve1}
    ${mont_curve1_point2} ${mont_curve1_point4})

do_test_curve_mul(elliptic_curve_montgomery_mul1 elliptic_mul_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} 3 ${mont_curve1_point3})
do_test_curve_mul(elliptic_curve_montgomery_mul2 elliptic_mul_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} 5 ${mont_curve1_point5})
do_test_curve_mul(elliptic_curve_montgomery_mul3 elliptic_mul_montgomery_affine ${mont_curve1}
    ${mont_curve1_point2} 3 ${mont_curve1_point6})
