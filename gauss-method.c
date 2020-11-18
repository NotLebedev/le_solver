#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss-method.h"
#include "debug.h"

/**
 * Применяет метод Гаусса к матрицам, результатом является "перемешанная треугольная матрица",
 * т.е. матрица которая отличается от треугольной только порядком столбоцв. Так как для всех
 * вычислений, кроме вычисления обратной матрицы порядок столбцов несущественен, то столбцы
 * матрицы остаются на местах, но их реальный порядок сохраняется дополнительно
 * @param a матрица системы
 * @param f столбец свободных членов (может быть NULL)
 * @param use_pivot выбирать наибольший по модулю элемент или произвольный (первый ненулевой)
 * элемент в качестве ведущего
 * @param operand дополнительная матрица над которой будут производиться те же действия что и над
 * a и f, используется для вычисления обратной матрицы (может быть NULL)
 * @return массив содержащий данные о перестановке столбцов или NULL в случае ошибки
 */
static size_t *elimination(Matrix *a, Matrix *f, _Bool use_pivot, Matrix *operand);
/**
 * Вычислить определитель приведённой к перемешанному треугольному виду матрицы
 * @param a матрица, должна быть приведена к перемешанному треугольному виду функцией elimination
 * @param col_order
 * @return
 */
static data_t determinant(Matrix *a, size_t *col_order);

/**
 * Проверить равно ли данное вещественное число нулю (с некоторой точностью)
 * @param d число для проверки
 * @return 1 если равно нулю, 0 иначе
 */
int eq_zero(data_t d)
{
#define eps (1e-300)
    return fabs(d) < eps;
#undef eps
}

GaussSolutionContainer *
o_gauss_solve(Matrix *a, Matrix *f, GaussSolverOptions *options)
{
    if (a == NULL || options == NULL || (options->solve_system && f == NULL) || a->col != a->row
            || (f != NULL && a->row != f->row)) {
#ifdef ASSERT_INCORRECT_ARGUMENTS
        fprintf(stderr, "Incorrect arguments passed to function gauss_solve\n");
#endif
        return NULL;
    }

    GaussSolutionContainer *result = calloc(1, sizeof(*result));
    if (result == NULL) {
        return NULL;
    }
    size_t *col_order = elimination(a, f, options->use_pivotal, NULL);
    if (col_order == NULL) {
        free(result);
        return NULL;
    }

    if (options->calc_determinant) {
        result->determinant = determinant(a, col_order);
    }

    if (options->solve_system) {

    }

    return result;
}

static size_t max_row_element(Matrix *a, const _Bool *cols_eliminated, size_t row);
static size_t nonzero_row_element(Matrix *a, const _Bool *cols_eliminated, size_t row);

static size_t *elimination(Matrix *a, Matrix *f, _Bool use_pivot, Matrix *operand)
{
    size_t *col_order = calloc(a->col, sizeof(*col_order));
    if (col_order == NULL) {
        return NULL;
    }
    _Bool *cols_eliminated = calloc(a->col, sizeof(*cols_eliminated));
    if (cols_eliminated == NULL) {
        free(col_order);
        return NULL;
    }
    for (size_t i = 0; i < a->row ; i++) {
        size_t pivot = use_pivot ? max_row_element(a, cols_eliminated, i) :
                nonzero_row_element(a, cols_eliminated, i);
        col_order[i] = pivot;
        cols_eliminated[pivot] = 1;
        data_t pivot_value = get_element(a, i, pivot);
        if (eq_zero(pivot_value)) {
            continue; // A row of all exactly zeroes should be skipped
        }

        for (size_t j = i + 1; j < a->row; j++) {
            data_t coefficient = - get_element(a, j, pivot) / pivot_value;
            mul_sub_row(a, i, j, coefficient);
            if (f != NULL) {
                mul_sub_row(f, i, j, coefficient);
            }
            if (operand != NULL) {
                mul_sub_row(operand, i, j, coefficient);
            }
        }

#ifdef PRINT_MATRIX_ELIMINATION_ITERATIONS
        print_matrix(stderr, a);
        if (f != NULL) {
            fputs("\n", stderr);
            print_matrix(stderr, f);
        }
        fputs("--------------------------------\n", stderr);
#endif
    }

    free(cols_eliminated);
    return col_order;
}

static size_t max_row_element(Matrix *a, const _Bool *cols_eliminated, size_t row)
{
    data_t max;
    size_t max_idx = 0;
    _Bool has_max = 0;
    for (size_t i = 0; i < a->col; i++) {
        data_t element = get_element(a, row, i);
        if (!cols_eliminated[i] && (!has_max || fabs(element) > max)) {
            max = fabs(element);
            max_idx = i;
            has_max = 1;
        }
    }

    return max_idx;
}

static size_t nonzero_row_element(Matrix *a, const _Bool *cols_eliminated, size_t row)
{
    for (size_t i = 0; i < a->col; i++) {
        if (!cols_eliminated[i] && !eq_zero(get_element(a, row, i))) {
            return i;
        }
    }

    return 0;
}

static data_t determinant(Matrix *a, size_t *col_order)
{
    data_t det = 1.0;
    for (size_t i = 0; i < a->row; i++) {
        det *= get_element(a, i, col_order[i]);
        for (size_t j = 0; j < i; j++) { // Multiply by (-1)^(number of inversions)
            if (col_order[j] > col_order[i]) {
                det *= -1.0;
            }
        }
    }

    return det;
}