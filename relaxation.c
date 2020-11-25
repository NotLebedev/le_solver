#include <math.h>

#include "relaxation.h"
#include "error.h"

Matrix *relaxation(Matrix *a, Matrix *f, data_t omega, data_t precision, size_t *iters, int *status) {
    if (a == NULL || f == NULL || a->col != a->row || f->row != a->row) {
        *status = INCORRECT_ARGS;
        return NULL;
    }

    // Матрица A должна удовлетворять условиям теоремы Самарского, т.е. быть самосопряженной, для этого домножим
    // систему уравнений с обоих сторон на A^T получая систему A^T * A * x = A^T * f
    Matrix *tm = transpose(a);
    if (tm == NULL) {
        *status = ALLOC_FAILED;
        return NULL;
    }
    if ((a = matrix_mul(tm, a)) == NULL) {
        *status = ALLOC_FAILED;
        free_matrix(tm);
        return NULL;
    }
    if ((f = matrix_mul(tm, f)) == NULL) {
        *status = ALLOC_FAILED;
        free_matrix(a);
        free_matrix(tm);
        return NULL;
    }
    free_matrix(tm);

    // Используется только один вектор решения, т.к. не требуется одновременно хранить x_i для более чем одной итерации,
    // т.к. x_i вычисляются последовательно
    Matrix *solution_cur = new_matrix(a->row, 1);
    if (solution_cur == NULL) {
        *status = ALLOC_FAILED;
        free_matrix(a);
        free_matrix(f);
        return NULL;
    }

    *iters = 0;
    data_t distance_squared = 0.0;
    do {
        distance_squared = 0.0; // Аккумулируем расстояние между векторами двух итераций нахождения решения
        for (size_t i = 0; i < a->row; i++) {
            data_t sum = 0.0; // Расчёт двух сум преобразуется в расчёт одной т.к. вектор-столбец решений
            // solution_cur содержит решения из вектора solution_prev начиная с i-ой строки
            for (size_t j = 0; j < a->row; j++) {
                sum += get_element(a, i, j) * get_element(solution_cur, j, 0);
            }

            // Находим x^(k+1)_i - x^k_i = w / a_(ii) * (f_i - sum)
            data_t variation = omega * (get_element(f, i, 0) - sum) / get_element(a, i, i);
            // Также добавляем (x^(k+1)_i - x^k_i) ^ 2 к расстоянию между между векторами двух
            // итераций нахождения решения
            distance_squared += pow(variation, 2);
            //
            set_element(solution_cur, i, 0, get_element(solution_cur, i, 0) + variation);
        }

        iters++;
    } while (sqrt(distance_squared) > precision); // Оценкой точности служит |x^(k+1) - x^k| < eps

    free_matrix(a);
    free_matrix(f);
    return solution_cur;
}
