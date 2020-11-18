#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"
#include "algs.h"

int
main(void)
{
    srand (124142);

    Matrix *a = new_matrix(10, 10);
    for (size_t i = 0; i < a->row; i++) {
        for (size_t j = 0; j < a->col; j++) {
            a->values[j + a->col * i] = (double)rand() / RAND_MAX * 10.0 - 5.0;
        }
    }

    Matrix *f = new_matrix(10, 1);
    for (size_t i = 0; i < a->row; i++) {
        f->values[i] = (double)rand() / RAND_MAX * 10.0 - 5.0;
    }

    print_matrix(stdout, a);
    print_matrix(stdout, f);

    int status;
    //printf("det = %" PR_DATA_T "\n", calc_determinant(a, 0, &status));
    //printf("det = %" PR_DATA_T "\n", calc_determinant(a, 1, &status));
    Matrix *answ = gauss_solve(a, f, 1, &status);
    print_matrix(stdout, answ);

    free_matrix(answ);
    free_matrix(a);
    free_matrix(f);
    return 0;
}
