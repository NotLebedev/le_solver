#ifndef LE_SOLVER_GAUSS_METHOD_H
#define LE_SOLVER_GAUSS_METHOD_H

#include "matrix.h"

typedef struct {
    _Bool solve_system; /* Вычислить вектро-столбец решений */
    _Bool use_pivotal; /* Использовать главный элемент в прямом ходе метода Гаусса */
    _Bool calc_determinant; /* Вычислить определитель матрицы системы */
    _Bool calc_inverse; /* Вычислить обратную матрицу к матрице системы */
} GaussSolverOptions;

typedef struct {
    Matrix *solution_vector; /* N x 1 вектор-столбец решений */
    data_t determinant; /* Определитель матрицы системы */
    Matrix inverse; /* Обратная матрица к матрице системы */
} GaussSolutionContainer;

/**
 * Solves the system Ax=f and calculates additional parameters as specified in options
 * @param a coefficient matrix of system (matrix A with dimensions N x N), will be modified in
 * process
 * @param f constant terms vector of system (vector-column f with dimensions N x 1), parameter is
 * ignored if solution of system is not calculated, will be modified in process
 * @param options structure specifying calculations to be performed
 * @return GaussSolutionContainer. NULL if an error happened or dimensions of matrices are incorrect
 */
GaussSolutionContainer *gauss_solve(Matrix *a, Matrix *f, GaussSolverOptions *options);

#endif //LE_SOLVER_GAUSS_METHOD_H
