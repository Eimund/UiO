/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  05.11.2014
 *
 *  c11 compiler
 */

#include "Matrix.h"

#define FLOAT double

int main() {
    auto matrix = Matrix<MatrixType::Tridiagonal_m1_C_m1, FLOAT>(4);
    matrix.Diagonal(2);
    MatrixCout(matrix);
    ArrayCout(matrix.factor, matrix.n);
    ArrayCout(matrix.n);
    //int hei = 0;
    return 0;
}
