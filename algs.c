#include <math.h>
#include <stdlib.h>

#include "algs.h"
#include "error.h"

static size_t *elimination(Matrix *a, Matrix *f, _Bool use_pivot);

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
    if ((a = copy_matrix(a)) == NULL) { // Копия для преобразований
        *status = ALLOC_FAILED;
        return 0.0;
    }
    size_t *col_order = elimination(a, NULL, use_pivot);
    if (col_order == NULL) {
        *status = ALLOC_FAILED;
        free_matrix(a);
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
    free_matrix(a);
    return det;
}

Matrix *gauss_solve(Matrix *a, Matrix *f, _Bool use_pivot, int *status) {
    if (a == NULL  || a->col != a->row || f == NULL) {
        *status = INCORRECT_ARGS;
        return NULL;
    }
    if ((a = copy_matrix(a)) == NULL) { // Копии для преобразований
        *status = ALLOC_FAILED;
        return NULL;
    }
    if ((f = copy_matrix(f)) == NULL) {
        *status = ALLOC_FAILED;
        free_matrix(a);
        return NULL;
    }

    size_t *col_order = elimination(a, f, use_pivot); // Выполняем прямой ход метода Гаусса
    if (col_order == NULL) {
        *status = ALLOC_FAILED;
        free_matrix(a);
        free_matrix(f);
        return NULL;
    }
    Matrix *answ = new_matrix(f->row, 1);
    if (answ == NULL) {
        *status = ALLOC_FAILED;
        free(col_order);
        free_matrix(a);
        free_matrix(f);
        return NULL;
    }

    // Обратный ход метода Гаусса
    for (size_t i = a->row; i > 0; i--) {
        data_t acc = get_element(f, i - 1, 0); // Находим сумму для всех найденных неизвестных и свободного члена
        for (size_t j = i; j < a->col; j++) {
            size_t col = col_order[j];
            acc -= get_element(a, i - 1, col) * get_element(answ, col, 0);
        }
        // Вычисляем неизвестную и записываем в вектор-столбец ответа
        set_element(answ, col_order[i - 1], 0, acc / get_element(a, i - 1, col_order[i - 1]));
    }

    free_matrix(a);
    free_matrix(f);
    return answ;
}

Matrix *calc_inverse(Matrix *a, int *status) {
    if (a == NULL  || a->col != a->row) {
        *status = INCORRECT_ARGS;
        return NULL;
    }
    if ((a = copy_matrix(a)) == NULL) { // Копия для преобразований
        *status = ALLOC_FAILED;
        return NULL;
    }
    // Создаём присоединёную единичную матрицу
    Matrix *res = new_matrix(a->row, a->col);
    if (res == NULL) {
        *status = ALLOC_FAILED;
        free_matrix(a);
        return NULL;
    }
    for (size_t i = 0; i < a->row; i++) {
        set_element(res, i, i, 1.0);
    }

    // В силу того, что алгоритм прямого хода метода Гаусса в функции elimination использует фиктивную перестановку
    // столбцов он не слишком удобен для нахождения обратной матрицы и не даёт выигрыша в производительности,
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
    }

    // Обратный ход, вычитаем из всех предыдущих строк обеих матриц текущую, домноженную на необходимый коэффициент
    for (size_t i = a->row; i > 0; i--) {
        for (size_t j = i - 1; j > 0; j--) {
            mul_sub_row(res, i - 1, j - 1, -get_element(a, j - 1, i - 1));
            mul_sub_row(a, i - 1, j - 1, -get_element(a, j - 1, i - 1));
        }
    }

    *status = OK;
    free_matrix(a);
    return res;
}

data_t matrix_norm(Matrix *a);

data_t calc_condition_number(Matrix *a, int *status)
{
    // Вычисляем число обусловленности как ||A|| * ||A^-1||
    data_t norm = matrix_norm(a);
    Matrix *inverse = calc_inverse(a, status);
    if (*status != OK) { // Не удалось вычислить обратную матрицу
        return 0.0;
    }
    return norm * matrix_norm(inverse);
}

data_t matrix_norm(Matrix *a)
{
    data_t max = 0.0;
    // Матричная норма ||A||_1
    for (size_t i = 0; i < a->row; i++) {
        data_t sum = 0.0; // Находим сумму модулей элементов строки
        for (size_t j = 0; j < a->col; j++) {
            sum += fabs(get_element(a, i, j));
        }

        // Проверям на максимальность и перезаписываем
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
    size_t *col_order = calloc(a->col, sizeof(*col_order));// Чтобы исбежать дополнительных вычислений вместо
    // перестановки столбцов при выборе главного элемента или в случае 0 диагонального сохраняется
    // дополнительный массив с порядком перестановки столбцов
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
        col_order[i] = pivot; // Записывем номер стобца в массив с порядком перестановки
        cols_eliminated[pivot] = 1; // Помечаем, что столбец был занулён
        data_t pivot_value = get_element(a, i, pivot);
        if (eq_zero(pivot_value)) {
            continue; // Если в качестве опорного эемента был выбран 0, то вся строка состоит из нулей и пропускается
        }

        for (size_t j = i + 1; j < a->row; j++) { // Из всех последующих строк вычитаем данную
            data_t coefficient = - get_element(a, j, pivot) / pivot_value; //  с необходимым коэффициентом
            mul_sub_row(a, i, j, coefficient);
            if (f != NULL) { // Функция может использоваться для вычисления определителя и тогда столбец свободных
                mul_sub_row(f, i, j, coefficient); // членов будет отсутствовать
            }
        }

    }

    free(cols_eliminated);
    return col_order;
}

static size_t max_row_element(Matrix *a, const _Bool *cols_eliminated, size_t row)
{
    data_t max = -1e300;
    size_t max_idx = 0;
    _Bool has_max = 0;
    for (size_t i = 0; i < a->col; i++) {
        data_t element = get_element(a, row, i);
        if (!cols_eliminated[i] && (!has_max || fabs(element) > max)) {
            max = fabs(element); // Максимальный элемент выбирается по модулю
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