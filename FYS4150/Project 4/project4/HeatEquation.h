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
    private: typedef HeatEquation<T, D> THIS;
    public: Property<PropertyType::ReadOnly,THIS,PropertySet<THIS,T>,T> theta;
    public: HeatEquation(unsigned int n) :
        Matrix<MatrixType::Tridiagonal_m1_C_m1, T>(n),
        theta(PROPERTY(this->theta, 0, Delegate<THIS,T,T>(this, &THIS::Theta))) {
        //theta(0, Delegate<THIS,T,T>(this, &THIS::Theta)) {
        this->Diagonal(2);
    }
    private: inline T Theta(T theta) {
        return theta;
    }
};

#endif // HEATEQUATION_H
