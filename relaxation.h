#ifndef LE_SOLVER_RELAXATION_H
#define LE_SOLVER_RELAXATION_H

#include "matrix.h"

Matrix *relaxation(Matrix *a, Matrix *f, data_t omega, data_t precision, int *status);

#endif //LE_SOLVER_RELAXATION_H
