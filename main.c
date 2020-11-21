#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"
#include "algs.h"

int
main(void)
{
    srand (13);

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

    printf("Матрица системы :\n");
    print_matrix(stdout, a);
    printf("\nМатрица-столбец коэффициентов :\n");
    print_matrix(stdout, f);

    int status;
    data_t det = calc_determinant(a, 0, &status);
    data_t det_pivot = calc_determinant(a, 1, &status);
    data_t condition_number = calc_condition_number(a, &status);
    Matrix *solution = gauss_solve(a, f, 0, &status);
    Matrix *solution_pivot = gauss_solve(a, f, 1, &status);
    Matrix *inverse = calc_inverse(a, &status);

    printf("\nРешение системы найденное методом Гаусса :\n");
    print_matrix(stdout, solution);
    printf("\nРешение системы найденное методом Гаусса с выбором главного элемента :\n");
    print_matrix(stdout, solution_pivot);
    printf("\nОбратная матрица : \n");
    print_matrix(stdout, inverse);

    printf("\nОпределитель вычисленный без выбора главного элемента : %" PR_DATA_T "\n", det);
    printf("Оперделитель вычисленный с выбором главного элемента : %" PR_DATA_T "\n", det_pivot);
    printf("Число обусловленности : %" PR_DATA_T "\n", condition_number);

    free_matrix(a);
    free_matrix(f);
    free_matrix(solution_pivot);
    free_matrix(solution);
    free_matrix(inverse);
    return 0;
}
