add_library(elliptic
    elliptic_curve_weierstrass_affine.c)

target_include_directories(elliptic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(${elliptic_include_dir} ${CMAKE_CURRENT_SOURCE_DIR})
