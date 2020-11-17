#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"
#include "algs.h"

int
main(void)
{
    srand (124142);

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

    int status;
    printf("det = %" PR_DATA_T "\n", calc_determinant(a, 0, &status));
    printf("det = %" PR_DATA_T "\n", calc_determinant(a, 1, &status));

    free(a);
    return 0;
}
