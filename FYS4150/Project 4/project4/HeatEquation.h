/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  07.11.2014
 *
 *  c11 compiler
 */
#ifndef HEATEQUATION_H
#define HEATEQUATION_H

#include "Matrix.h"

template<typename T, unsigned int D> class HeatEquation : public Matrix<MatrixType::Tridiagonal_m1_C_m1, T> {
    public: HeatEquation(unsigned int n) :
        Matrix<MatrixType::Tridiagonal_m1_C_m1, T>(n) {
        this->Diagonal(2);
    }
};

#endif // HEATEQUATION_H
