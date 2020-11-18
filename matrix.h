#ifndef LE_SOLVER_MATRIX_H
#define LE_SOLVER_MATRIX_H

#include <stddef.h>
#include <stdio.h>

typedef double data_t;
#define PR_DATA_T "lf"

typedef struct
{
    size_t row;
    size_t col;
    data_t values[];
} Matrix;

Matrix *new_matrix(size_t row, size_t col);
void free_matrix(Matrix *matrix);

data_t get_element(Matrix *matrix, size_t row, size_t col);
void set_element(Matrix *matrix, size_t row, size_t col, data_t val);
/**
 * Прибавить к строке d матрицы a строку s, умноженную на c
 * @param matrix матрица для преобразования
 */
void mul_sub_row(Matrix *matrix, size_t s, size_t d, data_t c);

void mul_row(Matrix *matrix, size_t row, data_t c);

void swap_row(Matrix *matrix, size_t s, size_t d);

void print_matrix(FILE *file, Matrix *matrix);

#endif //LE_SOLVER_MATRIX_H
