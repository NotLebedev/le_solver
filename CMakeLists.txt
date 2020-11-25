cmake_minimum_required(VERSION 3.16)
project(le_solver C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -std=gnu11 -fsanitize=undefined -Wall -Werror -Wno-pointer-sign -Wformat -Wformat-overflow -Wformat-security -Wnull-dereference -Wignored-qualifiers -Wshift-negative-value -Wswitch-default -Wduplicated-branches -Wduplicated-branches -Wfloat-equal -Wshadow -Wpointer-arith -Wpointer-compare -Wtype-limits -Wwrite-strings -Wdangling-else -Wempty-body -Wlogical-op -Wstrict-prototypes -Wold-style-declaration -Wold-style-definition -Wmissing-parameter-type -Wmissing-field-initializers -Wnested-externs -Wvla-larger-than=4096 -Wno-unused-result -lm")

add_executable(le_solver main.c matrix.h matrix.c algs.c algs.h error.h relaxation.c relaxation.h)
target_link_libraries(le_solver m)