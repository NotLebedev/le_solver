#ifndef LE_SOLVER_ALGS_H
#define LE_SOLVER_ALGS_H

#include "matrix.h"

Matrix *gauss_solve(Matrix *a, Matrix *f, _Bool use_pivot);

data_t calc_determinant(Matrix *a, _Bool use_pivot);

Matrix *calc_inverse(Matrix *a);

data_t calc_condition_number(Matrix *a);

#endif //LE_SOLVER_ALGS_H
