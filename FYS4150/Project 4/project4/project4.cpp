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
    diffusion.len1 = 20;
    diffusion.len1 = &diffusion.len2;
    //diffusion.n = diffusion.n;
    ArrayCout(&diffusion.n);
    int hei = 0;
    return 0;
}
