#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "matrix.h"
#include "debug.h"

Matrix *
new_matrix(size_t row, size_t col)
{
    Matrix *matrix = calloc(1, sizeof(*matrix) + sizeof(data_t) * row * col);
    if (matrix == NULL) {
        return NULL;
    }
    matrix->col = col;
    matrix->row = row;
    return matrix;
}

void
free_matrix(Matrix *matrix)
{
    free(matrix);
}

Matrix *copy_matrix(Matrix *matrix) {
    Matrix *copy = new_matrix(matrix->row, matrix->col);
    *copy = *matrix;
    memcpy(copy->values, matrix->values, matrix->col * matrix->row * sizeof(*(copy->values)));
    return copy;
}

data_t
get_element(Matrix *matrix, size_t row, size_t col)
{
#ifdef ASSERT_MATRIX_RANGES
    assert(row < matrix->row);
    assert(col < matrix->col);
#endif
    return matrix->values[col + matrix->col * row];
}

void
set_element(Matrix *matrix, size_t row, size_t col, data_t val)
{
#ifdef ASSERT_MATRIX_RANGES
    assert(row < matrix->row);
    assert(col < matrix->col);
#endif
    matrix->values[col + matrix->col * row] = val;
}

void
mul_sub_row(Matrix *matrix, size_t s, size_t d, data_t c)
{
#ifdef ASSERT_MATRIX_RANGES
    assert(s < matrix->row);
    assert(d < matrix->row);
#endif
    for (size_t i = 0; i < matrix->col; i++) {
        matrix->values[i + d * matrix->col] += c * matrix->values[i + s * matrix->col];
    }
}

void
print_matrix(FILE *file, Matrix *matrix)
{
#ifdef PRINT_MATRIX_HUMAN_READABLE
    for (size_t i = 0; i < matrix->row; i++) {
        for (size_t j = 0; j < matrix->col; j++) {
            fprintf(file, "%10.5" PR_DATA_T " ", matrix->values[j + matrix->col * i]);
        }
        fputs("\n", file);
    }
#else
    fputs("{", file);
    for (size_t i = 0; i < matrix->row; i++) {
        fputs("{", file);
        for (size_t j = 0; j < matrix->col; j++) {
            fprintf(file, "%" PR_DATA_T "%s", matrix->values[j + matrix->col * i],
                    j != matrix-> col - 1 ? ", " : "");
        }
        fputs(i != matrix->row - 1 ? "}," : "}", file);
    }
    fputs("}\n", file);
#endif
}

void swap_row(Matrix *matrix, size_t s, size_t d) {
#ifdef ASSERT_MATRIX_RANGES
    assert(s < matrix->row);
    assert(d < matrix->row);
#endif
    for (size_t i = 0; i < matrix->col; i++) {
        data_t tmp = matrix->values[i + d * matrix->col];
        matrix->values[i + d * matrix->col] = matrix->values[i + s * matrix->col];
        matrix->values[i + s * matrix->col] = tmp;
    }
}

void mul_row(Matrix *matrix, size_t row, data_t c) {
#ifdef ASSERT_MATRIX_RANGES
    assert(row < matrix->row);
#endif
    for (size_t i = 0; i < matrix->col; i++) {
        matrix->values[i + row * matrix->col] *= c;
    }
}