#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"
#include "gauss-method.h"

int
main(void)
{
    srand ( 12123);

    Matrix *a = new_matrix(5, 5);
    for (size_t i = 0; i < a->row; i++) {
        for (size_t j = 0; j < a->col; j++) {
            a->values[j + a->col * i] = (double)rand() / RAND_MAX * 10.0 - 5.0;
        }
    }

    Matrix *f = new_matrix(5, 1);
    for (size_t i = 0; i < a->row; i++) {
        f->values[i] = (double)rand() / RAND_MAX * 10.0 - 5.0;
    }

    print_matrix(stdout, a);

    GaussSolverOptions opt = {0};
    opt.use_pivotal = 1;
    opt.calc_determinant = 1;
    GaussSolutionContainer *res = gauss_solve(a, f, &opt);

    printf("det = %" PR_DATA_T "\n", res->determinant);

    free(res);
    free(a);
    return 0;
}
