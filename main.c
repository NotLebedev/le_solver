#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "algs.h"
#include "relaxation.h"

extern int mode;

data_t gen_function_a(size_t n, size_t m, size_t i, size_t j);
data_t gen_function_f(size_t n, size_t m, size_t i);

/*
 * Первый аргумент -- число 1 или 2 для выбора между методом Гаусса и методом верхней релакскации, соответственно
 * Второй аргумень -- число 1 или 2 для выбора формата вывода матрицы -- человекочитаемо или машинописно, соответственно
 * Третий аргумент -- число 1 или 2 для выбора ввода матрицы вручную или генерации при помощи функции
 * Для случая задания матрицы вручную передаётся число n -- размер матрицы A и вектора-столбца f, которые
 * задаются на стандартном потоке ввода
 * Для случая генерации матрицы при помощи функции задаются два числа n и m
 */
int main(int argc, char *argv[]) {
    size_t n = strtoul(argv[4], NULL, 0);
    mode = strtoul(argv[2], NULL, 0) - 1;

    Matrix *a = new_matrix(n, n);
    if (a == NULL) {
        fprintf(stderr, "Ошибка выделения памяти\n");
        return 1;
    }
    Matrix *f = new_matrix(n, 1);
    if (f == NULL) {
        fprintf(stderr, "Ошибка выделения памяти\n");
        free_matrix(a);
        return 1;
    }

    if (strtoul(argv[3], NULL, 0) == 1) {

    } else if (strtoul(argv[3], NULL, 0) == 2) {
        size_t m = strtoul(argv[5], NULL, 0);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                set_element(a, i, j, gen_function_a(n, m, i, j));
            }
        }
        for (size_t i = 0; i < n; i++) {
            set_element(f, i, 0, gen_function_f(n, m, i));
        }
    } else {
        fprintf(stderr, "Некорректный режим задания матрицы\n");
        free_matrix(a);
        free_matrix(f);
        return 1;
    }

    printf("Матрица системы :\n");
    print_matrix(stdout, a);
    printf("\nМатрица-столбец коэффициентов :\n");
    print_matrix(stdout, f);

    if (strtoul(argv[1], NULL, 0) == 1) {
        int status;
        data_t det = calc_determinant(a, 0, &status);
        data_t det_pivot = calc_determinant(a, 1, &status);
        data_t condition_number = calc_condition_number(a, &status);
        Matrix *solution = gauss_solve(a, f, 0, &status);
        Matrix *solution_pivot = gauss_solve(a, f, 1, &status);
        Matrix *inverse = calc_inverse(a, &status);
        if (status < 0) {
            fprintf(stderr, "Ошибка во время выполнения программы\n");
            return 1;
        }

        printf("\nРешение системы найденное методом Гаусса :\n");
        print_matrix(stdout, solution);
        printf("\nРешение системы найденное методом Гаусса с выбором главного элемента :\n");
        print_matrix(stdout, solution_pivot);
        printf("\nОбратная матрица : \n");
        print_matrix(stdout, inverse);

        printf("\nОпределитель вычисленный без выбора главного элемента : %" PR_DATA_T "\n", det);
        printf("Оперделитель вычисленный с выбором главного элемента : %" PR_DATA_T "\n", det_pivot);
        printf("Число обусловленности : %" PR_DATA_T "\n", condition_number);

        free_matrix(solution_pivot);
        free_matrix(solution);
        free_matrix(inverse);
    } else if (strtoul(argv[1], NULL, 0) == 2) {
        int status;
        size_t iter;
        print_matrix(stdout, relaxation(a, f, 1, 1e-10, &iter, &status));
    }

    free_matrix(a);
    free_matrix(f);
    return 0;
}

data_t gen_function_a(size_t n, size_t m, size_t i, size_t j)
{
    return i == j ? n + ((data_t) m) * m + ((data_t) j) / m + ((data_t) i) / n
            : ((data_t) i + j) / (((data_t)m) + n);
}

data_t gen_function_f(size_t n, size_t m, size_t i)
{
    return ((data_t) m) * i + n;
}