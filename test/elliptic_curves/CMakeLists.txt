macro(prepare_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} elliptic m gmp)
    target_include_directories(${target} PUBLIC ${elliptic_include_dir})
endmacro(prepare_executable)

### Macros for creating tests with affine coordinates
macro(do_test_curve_sum_affine test_name test_executable
        c_a c_b c_m p1_x p1_y p2_x p2_y r_x r_y)
    add_test (${test_name} ${test_executable}
        ${c_a} ${c_b} ${c_m} ${p1_x} ${p1_y} ${p2_x} ${p2_y})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION "${r_x}\t${r_y}")
endmacro(do_test_curve_sum_affine)

macro(do_test_curve_double_affine test_name test_executable
        c_a c_b c_m p1_x p1_y r_x r_y)
    add_test (${test_name} ${test_executable}
        ${c_a} ${c_b} ${c_m} ${p1_x} ${p1_y})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION "${r_x}\t${r_y}")
endmacro(do_test_curve_double_affine)

macro(do_test_curve_mul_affine test_name test_executable
        c_a c_b c_m p1_x p1_y times r_x r_y)
    add_test (${test_name} ${test_executable}
        ${c_a} ${c_b} ${c_m} ${p1_x} ${p1_y} ${times})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION "${r_x}\t${r_y}")
endmacro(do_test_curve_mul_affine)


### Macros for creating tests with projective coordinates
# There are many points with different coordinates that are
# congruent:
# x1/z1 mod m === x2/z2 mod m
# so, all the points returned by the tests should have
# z = 1, so the x and y coordinates are the same
macro(do_test_curve_sum_projective test_name test_executable
        c_a c_b c_m p1_x p1_y p1_z p2_x p2_y p2_z r_x r_y)
    add_test (${test_name} ${test_executable}
        ${c_a} ${c_b} ${c_m} ${p1_x} ${p1_y} ${p1_z} ${p2_x} ${p2_y} ${p2_z})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION "${r_x}\t${r_y}")
endmacro(do_test_curve_sum_projective)

macro(do_test_curve_double_projective test_name test_executable
        c_a c_b c_m p1_x p1_y p1_z r_x r_y)
    add_test (${test_name} ${test_executable}
        ${c_a} ${c_b} ${c_m} ${p1_x} ${p1_y} ${p1_z})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION "${r_x}\t${r_y}")
endmacro(do_test_curve_double_projective)

macro(do_test_curve_mul_projective test_name test_executable
        c_a c_b c_m p1_x p1_y p1_z times r_x r_y)
    add_test (${test_name} ${test_executable}
        ${c_a} ${c_b} ${c_m} ${p1_x} ${p1_y} ${p1_z} ${times})
    set_tests_properties (${test_name}
        PROPERTIES PASS_REGULAR_EXPRESSION "${r_x}\t${r_y}")
endmacro(do_test_curve_mul_projective)


prepare_executable(elliptic_sum_weierstrass_affine    elliptic_sum_weierstrass_affine.c)
prepare_executable(elliptic_double_weierstrass_affine elliptic_double_weierstrass_affine.c)
prepare_executable(elliptic_mul_weierstrass_affine    elliptic_mul_weierstrass_affine.c)

prepare_executable(elliptic_sum_weierstrass_projective    elliptic_sum_weierstrass_projective.c)
prepare_executable(elliptic_double_weierstrass_projective elliptic_double_weierstrass_projective.c)
prepare_executable(elliptic_mul_weierstrass_projective    elliptic_mul_weierstrass_projective.c)

prepare_executable(elliptic_sum_montgomery_affine    elliptic_sum_montgomery_affine.c)
prepare_executable(elliptic_double_montgomery_affine elliptic_double_montgomery_affine.c)
prepare_executable(elliptic_mul_montgomery_affine    elliptic_mul_montgomery_affine.c)


### Tests for Weierstrass curves with affine coordinates
set(weir_curve1 5 -5 65537)
set(weir_curve1_point1     1     1)
set(weir_curve1_point2    14 65484) # 2 * (1,1)
set(weir_curve1_point3 24821 62767) # 3 * (1,1)
set(weir_curve1_point4 42740 22681) # 4 * (1,1)
set(weir_curve1_point5 50687 11143) # 5 * (1,1)
set(weir_curve1_point6 45879 59979) # 6 * (1,1)

do_test_curve_sum_affine(elliptic_curve_weierstrass_affine_sum1 elliptic_sum_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} ${weir_curve1_point1} ${weir_curve1_point2})
do_test_curve_sum_affine(elliptic_curve_weierstrass_affine_sum2 elliptic_sum_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} ${weir_curve1_point2} ${weir_curve1_point3})
do_test_curve_sum_affine(elliptic_curve_weierstrass_affine_sum3 elliptic_sum_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} ${weir_curve1_point3} ${weir_curve1_point4})
do_test_curve_sum_affine(elliptic_curve_weierstrass_affine_sum4 elliptic_sum_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point2} ${weir_curve1_point3} ${weir_curve1_point5})

do_test_curve_double_affine(elliptic_curve_weierstrass_affine_double1 elliptic_double_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} ${weir_curve1_point2})
do_test_curve_double_affine(elliptic_curve_weierstrass_affine_double2 elliptic_double_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point2} ${weir_curve1_point4})

do_test_curve_mul_affine(elliptic_curve_weierstrass_affine_mul1 elliptic_mul_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} 3 ${weir_curve1_point3})
do_test_curve_mul_affine(elliptic_curve_weierstrass_affine_mul2 elliptic_mul_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point1} 5 ${weir_curve1_point5})
do_test_curve_mul_affine(elliptic_curve_weierstrass_affine_mul3 elliptic_mul_weierstrass_affine
    ${weir_curve1} ${weir_curve1_point2} 3 ${weir_curve1_point6})


### Tests for Weierstrass curves with projective coordinates
set(weir_curve2 -111 -110 65537)

set(weir_curve2_point1    15    40     1)
set(weir_curve2_point2 60519 33020 53241) # 2 * (15,40,1)
set(weir_curve2_point3 37252 15872 46509) # 3 * (15,40,1)
set(weir_curve2_point4 38881  4337 45041) # 4 * (15,40,1)

set(weir_curve2_point2_normal 14274 53446) # 2 * (15,40,1)
set(weir_curve2_point3_normal 15690 43989) # 3 * (15,40,1)
set(weir_curve2_point4_normal 39037 26095) # 4 * (15,40,1)
set(weir_curve2_point5_normal 61487 23562) # 4 * (15,40,1)
set(weir_curve2_point6_normal  2366 21609) # 5 * (15,40,1)

do_test_curve_sum_projective(
    elliptic_curve_weierstrass_projective_sum1
    elliptic_sum_weierstrass_projective
    ${weir_curve2} ${weir_curve2_point1} ${weir_curve2_point1} ${weir_curve2_point2_normal})
do_test_curve_sum_projective(
    elliptic_curve_weierstrass_projective_sum2
    elliptic_sum_weierstrass_projective
    ${weir_curve2} ${weir_curve2_point1} ${weir_curve2_point2} ${weir_curve2_point3_normal})
do_test_curve_sum_projective(
    elliptic_curve_weierstrass_projective_sum3
    elliptic_sum_weierstrass_projective
    ${weir_curve2} ${weir_curve2_point1} ${weir_curve2_point3} ${weir_curve2_point4_normal})
do_test_curve_sum_projective(
    elliptic_curve_weierstrass_projective_sum4
    elliptic_sum_weierstrass_projective
    ${weir_curve2} ${weir_curve2_point2} ${weir_curve2_point3} ${weir_curve2_point5_normal})

do_test_curve_double_projective(
    elliptic_curve_weierstrass_projective_double1
    elliptic_double_weierstrass_projective
    ${weir_curve2} ${weir_curve2_point1} ${weir_curve2_point2_normal})
do_test_curve_double_projective(
    elliptic_curve_weierstrass_projective_double2
    elliptic_double_weierstrass_projective
    ${weir_curve2} ${weir_curve2_point2} ${weir_curve2_point4_normal})

do_test_curve_mul_projective(
    elliptic_curve_weierstrass_projective_mul1
    elliptic_mul_weierstrass_projective
    ${weir_curve2} ${weir_curve2_point1} 3 ${weir_curve2_point3_normal})
do_test_curve_mul_projective(
    elliptic_curve_weierstrass_projective_mul2
    elliptic_mul_weierstrass_projective
    ${weir_curve2} ${weir_curve2_point1} 5 ${weir_curve2_point5_normal})
do_test_curve_mul_projective(
    elliptic_curve_weierstrass_projective_mul3
    elliptic_mul_weierstrass_projective
    ${weir_curve2} ${weir_curve2_point2} 3 ${weir_curve2_point6_normal})


### Tests for Montgomery curves with affine coordinates
set(mont_curve1 5 3 65537)
set(mont_curve1_point1      3     5)
set(mont_curve1_point2   6117  1282) # 2 * (3,5)
set(mont_curve1_point3  54943 58216) # 3 * (3,5)
set(mont_curve1_point4  62979 59844) # 4 * (3,5)
set(mont_curve1_point5  23228 18695) # 5 * (3,5)
set(mont_curve1_point6  12508 22665) # 6 * (3,5)

do_test_curve_sum_affine(elliptic_curve_montgomery_affine_sum1 elliptic_sum_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} ${mont_curve1_point1} ${mont_curve1_point2})
do_test_curve_sum_affine(elliptic_curve_montgomery_affine_sum2 elliptic_sum_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} ${mont_curve1_point2} ${mont_curve1_point3})
do_test_curve_sum_affine(elliptic_curve_montgomery_affine_sum3 elliptic_sum_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} ${mont_curve1_point3} ${mont_curve1_point4})
do_test_curve_sum_affine(elliptic_curve_montgomery_affine_sum4 elliptic_sum_montgomery_affine ${mont_curve1}
    ${mont_curve1_point2} ${mont_curve1_point3} ${mont_curve1_point5})

do_test_curve_double_affine(elliptic_curve_montgomery_affine_double1 elliptic_double_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} ${mont_curve1_point2})
do_test_curve_double_affine(elliptic_curve_montgomery_affine_double2 elliptic_double_montgomery_affine ${mont_curve1}
    ${mont_curve1_point2} ${mont_curve1_point4})

do_test_curve_mul_affine(elliptic_curve_montgomery_affine_mul1 elliptic_mul_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} 3 ${mont_curve1_point3})
do_test_curve_mul_affine(elliptic_curve_montgomery_affine_mul2 elliptic_mul_montgomery_affine ${mont_curve1}
    ${mont_curve1_point1} 5 ${mont_curve1_point5})
do_test_curve_mul_affine(elliptic_curve_montgomery_affine_mul3 elliptic_mul_montgomery_affine ${mont_curve1}
    ${mont_curve1_point2} 3 ${mont_curve1_point6})
