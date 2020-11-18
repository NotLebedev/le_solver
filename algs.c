#include <math.h>
#include <stdlib.h>

#include "algs.h"
#include "error.h"
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
static size_t *elimination(Matrix *a, Matrix *f, _Bool use_pivot);

/**
* Проверить равно ли данное вещественное число нулю (с некоторой точностью)
* @param d число для проверки
* @return 1 если равно нулю, 0 иначе
*/
static int eq_zero(data_t d)
{
#define eps (1e-300)
    return fabs(d) < eps;
#undef eps
}

data_t calc_determinant(Matrix *a, _Bool use_pivot, int *status) {
    if (a == NULL  || a->col != a->row) {
        *status = INCORRECT_ARGS;
        return 0.0;
    }

    size_t *col_order = elimination(a, NULL, use_pivot);
    if (col_order == NULL) {
        *status = ALLOC_FAILED;
        return 0.0;
    }

    data_t det = 1.0;
    for (size_t i = 0; i < a->row; i++) {
        det *= get_element(a, i, col_order[i]); // Перемножаем диагональные с точности до перестановки столбцов элементы
        for (size_t j = 0; j < i; j++) { // Домножаем на (-1)^(количество инверсий)
            if (col_order[j] > col_order[i]) {
                det *= -1.0;
            }
        }
    }

    *status = OK;
    free(col_order);
    return det;
}

Matrix *gauss_solve(Matrix *a, Matrix *f, _Bool use_pivot, int *status) {
    if (a == NULL  || a->col != a->row || f == NULL) {
        *status = INCORRECT_ARGS;
        return NULL;
    }

    size_t *col_order = elimination(a, f, use_pivot);
    if (col_order == NULL) {
        *status = ALLOC_FAILED;
        return NULL;
    }

    Matrix *answ = new_matrix(f->row, 1);
    if (answ == NULL) {
        *status = ALLOC_FAILED;
        free(col_order);
        return NULL;
    }

    // Обратный ход метода Гаусса
    for (size_t i = a->row; i > 0; i--) {
        data_t acc = get_element(f, i - 1, 0); // Находим сумму для всех найденных неизвестных и свободного члена
        for (size_t j = i; j < a->col; j++) {
            size_t col = col_order[j];
            acc -= get_element(a, i - 1, col) * get_element(answ, col, 0);
        }

        set_element(answ, col_order[i - 1], 0, acc / get_element(a, i - 1, col_order[i - 1]));
    }

    return answ;
}

Matrix *calc_inverse(Matrix *a, int *status) {
    if (a == NULL  || a->col != a->row) {
        *status = INCORRECT_ARGS;
        return NULL;
    }

    // Создаём присоединёную единичную матрицу
    Matrix *res = new_matrix(a->row, a->col);
    if (res == NULL) {
        *status = ALLOC_FAILED;
        return NULL;
    }
    for (size_t i = 0; i < a->row; i++) {
        set_element(res, i, i, 1.0);
    }

    // В силу того, что алгоритм прямого хода метода Гаусса в функции elemination использует фиктивную перестановку
    // столбцов он не слишком удобен для нахождения обратной матрицы и не даёт выигрыша в производительности
    // потому используется другая реализация прямого хода
    for (size_t i = 0; i < a->row; i++) {
        if (eq_zero(get_element(a, i, i))) { // Диагональный элемент необходимо выбрать ненулевым
            for (size_t j = i + 1; j < a->row; j++) {
                if (!eq_zero(get_element(a, j, i))) { // Находим строку с ненулевым элементом
                    swap_row(a, i, j); // Меняем текущую строку с найденной в обеих матрицах
                    swap_row(res, i, j);
                }
            }
        }
        mul_row(res, i, 1.0 / get_element(a, i, i)); // Делим строку обеих матриц на диагональный элемент
        mul_row(a, i, 1.0 / get_element(a, i, i));
        for (size_t j = i + 1; j < a->row; j++) {
            mul_sub_row(res, i, j, -get_element(a, j, i)); // Вычитаем из всех последующих строк обеих матриц данную
            mul_sub_row(a, i, j, -get_element(a, j, i)); // домноженную на необходимый коэффициент
        }
#ifdef PRINT_MATRIX_ELIMINATION_ITERATIONS
        print_matrix(stderr, a);
        fputs("\n", stderr);
        print_matrix(stderr, res);
        fputs("--------------------------------\n", stderr);
#endif
    }

    // Обратный ход, вычитаем из всех предыдущих строк обеих матриц текущую, домноженную на необходимый коэффициент
    for (size_t i = a->row; i > 0; i--) {
        for (size_t j = i - 1; j > 0; j--) {
            mul_sub_row(res, i - 1, j - 1, -get_element(a, j - 1, i - 1));
            mul_sub_row(a, i - 1, j - 1, -get_element(a, j - 1, i - 1));
        }
    }

    *status = OK;
    return res;
}

data_t matrix_norm(Matrix *a);

data_t calc_condition_number(Matrix *a, int *status)
{
    data_t norm = matrix_norm(a);
    fprintf(stderr, "Norm(A) = %" PR_DATA_T "\n", norm);
    Matrix *inverse = calc_inverse(a, status);
    print_matrix(stderr, inverse);
    if (*status != OK) {
        return 0.0;
    }
    fprintf(stderr, "Norm(A) = %" PR_DATA_T "\n", matrix_norm(inverse));
    return norm * matrix_norm(inverse);
}

data_t matrix_norm(Matrix *a)
{
    data_t max = 0.0;
    for (size_t i = 0; i < a->row; i++) {
        data_t sum = 0.0;
        for (size_t j = 0; j < a->col; j++) {
            sum += fabs(get_element(a, i, j));
        }

        if (sum > max) {
            max = sum;
        }
    }
    return max;
}

static size_t max_row_element(Matrix *a, const _Bool *cols_eliminated, size_t row);
static size_t nonzero_row_element(Matrix *a, const _Bool *cols_eliminated, size_t row);

static size_t *elimination(Matrix *a, Matrix *f, _Bool use_pivot)
{
    size_t *col_order = calloc(a->col, sizeof(*col_order));// Чтобы исбежать дополнительных вычислений вместо перестановки
    // столбцов при выборе главного элемента или в случае 0 диагонального сохраняется дополнительный массив с порядком
    // перестановки столбцов
    if (col_order == NULL) {
        return NULL;
    }
    _Bool *cols_eliminated = calloc(a->col, sizeof(*cols_eliminated)); // Также сохраняем информацию о том, какие
    // столбцы были уже занулены каждом шаге
    if (cols_eliminated == NULL) {
        free(col_order);
        return NULL;
    }
    for (size_t i = 0; i < a->row ; i++) {
        size_t pivot = use_pivot ? max_row_element(a, cols_eliminated, i) : // В зависимости от режима работы
                       nonzero_row_element(a, cols_eliminated, i); // выбирается либо максимальный элемент либо первый
                       // ненулевой в строке
        col_order[i] = pivot;
        cols_eliminated[pivot] = 1;
        data_t pivot_value = get_element(a, i, pivot);
        if (eq_zero(pivot_value)) {
            continue; // Если в качестве опорного эемента был выбран 0, то вся строка состоит из нулей и пропускается
        }

        for (size_t j = i + 1; j < a->row; j++) {
            data_t coefficient = - get_element(a, j, pivot) / pivot_value;
            mul_sub_row(a, i, j, coefficient);
            if (f != NULL) { // Функция может использоваться для вычисления определителя и тогда столбец свободных
                mul_sub_row(f, i, j, coefficient); // членов будет отсутствовать
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