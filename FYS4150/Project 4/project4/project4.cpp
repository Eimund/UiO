/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  05.11.2014
 *
 *  c11 compiler
 */

#include "HeatEquation.h"

#define FLOAT double

int main() {
    auto diffusion = HeatEquation<FLOAT, 1>(0.2, ARRAYLIST(FLOAT,0.9,0.1), ARRAYLIST(FLOAT,9.76,1.45), ARRAYLIST(int,4,6));
    diffusion.lower[0] = 0.8;
    diffusion.upper[0] = 8.1;
    diffusion.lower[1] = 0.23;
    diffusion.upper[1] = 2.13;
    MatrixCout(*diffusion.matrix[0]);
    diffusion.u[0][0] = 1;
    ArrayCout(& &diffusion.X[1]);
    ArrayCout(diffusion.u[0], diffusion.n[1]);
    //auto hei = 0;
    return 0;
}
