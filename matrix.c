#include <stdlib.h>
#include <string.h>

#include "matrix.h"

Matrix *new_matrix(size_t row, size_t col) {
    Matrix *matrix = calloc(1, sizeof(*matrix) + sizeof(data_t) * row * col);
    if (matrix == NULL) {
        return NULL;
    }
    matrix->col = col;
    matrix->row = row;
    return matrix;
}

void free_matrix(Matrix *matrix) {
    free(matrix);
}

Matrix *copy_matrix(Matrix *matrix) {
    Matrix *copy = new_matrix(matrix->row, matrix->col);
    if (copy == NULL) {
        return NULL;
    }
    *copy = *matrix;
    memcpy(copy->values, matrix->values, matrix->col * matrix->row * sizeof(*(copy->values)));
    return copy;
}

data_t get_element(Matrix *matrix, size_t row, size_t col) {
    return matrix->values[col + matrix->col * row];
}

void set_element(Matrix *matrix, size_t row, size_t col, data_t val) {
    matrix->values[col + matrix->col * row] = val;
}

void mul_sub_row(Matrix *matrix, size_t s, size_t d, data_t c) {
    for (size_t i = 0; i < matrix->col; i++) {
        matrix->values[i + d * matrix->col] += c * matrix->values[i + s * matrix->col];
    }
}

int mode = 0;

void print_matrix(FILE *file, Matrix *matrix) {
    if (mode == 0) {
        for (size_t i = 0; i < matrix->row; i++) {
            for (size_t j = 0; j < matrix->col; j++) {
                fprintf(file, "%10.5" PR_DATA_T " ", matrix->values[j + matrix->col * i]);
            }
            fputs("\n", file);
        }
    } else {
        fputs("{", file);
        for (size_t i = 0; i < matrix->row; i++) {
            fputs("{", file);
            for (size_t j = 0; j < matrix->col; j++) {
                fprintf(file, "%" PR_DATA_T "%s", matrix->values[j + matrix->col * i],
                        j != matrix->col - 1 ? ", " : "");
            }
            fputs(i != matrix->row - 1 ? "}," : "}", file);
        }
        fputs("}\n", file);
    }
}

void swap_row(Matrix *matrix, size_t s, size_t d) {
    for (size_t i = 0; i < matrix->col; i++) {
        data_t tmp = matrix->values[i + d * matrix->col];
        matrix->values[i + d * matrix->col] = matrix->values[i + s * matrix->col];
        matrix->values[i + s * matrix->col] = tmp;
    }
}

void mul_row(Matrix *matrix, size_t row, data_t c) {
    for (size_t i = 0; i < matrix->col; i++) {
        matrix->values[i + row * matrix->col] *= c;
    }
}

Matrix *matrix_mul(Matrix *a, Matrix *b) {
    Matrix *res = new_matrix(a->row, b->col);
    if (res == NULL) {
        return NULL;
    }
    for (size_t i = 0; i < a->row; i++) {
        for (size_t j = 0; j < b->col; j++) {
            for (size_t k = 0; k < a->col; k++) {
                res->values[j + res->col * i] += a->values[k + a->col * i] * b->values[j + b->col * k];
            }
        }
    }
    return res;
}

Matrix *transpose(Matrix *matrix) {
    Matrix *res = new_matrix(matrix->col, matrix->row);
    if (res == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < res->row; i++) {
        for (size_t j = 0; j < res->col; j++) {
            res->values[j + res->col * i] = matrix->values[i + matrix->col * j];
        }
    }
    return res;
}

