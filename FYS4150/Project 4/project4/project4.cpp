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
    auto diffusion = HeatEquation<FLOAT, 1>(4);
    MatrixCout(diffusion);
    diffusion.n = 10;
    diffusion.theta = 0.33;
    cout << diffusion.theta << endl;
    ArrayCout(&diffusion.n);
    int hei = 0;
    return 0;
}
