cmake_minimum_required(VERSION 2.8)
project(Math_stuff)

set(CMAKE_BUILD_TYPE Debug)

# Only set up to -O2 in release mode as -O3 may break
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG"
    CACHE STRING "Flags used by the compiler during release builds." FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG"
    CACHE STRING "Flags used by the compiler during release builds." FORCE)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra -pedantic")

add_subdirectory(src)
