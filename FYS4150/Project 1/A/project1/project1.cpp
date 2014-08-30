/*
 *  FYS4150 - Computational Physics - Project 1
 *
 *  Written by: Eimund Smestad
 *
 *  30.08.2014
 */

#include <cmath>
#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;

#define ARRAY_SIZE(X) sizeof(X)/sizeof(*X)      // Calculate static array size at compile time

template<class T> T* CopyOfArray(T* arr, unsigned int n);
template<class T> T* TridiagonalSolutionGeneral(T* a, T*b, T* c, T* f, unsigned int n);

int main() {
    unsigned int n[] = {4,100,1000};           // Size of matrix to solve
    double _x = 0, x_ = 1;                      // Solution interval

    for(unsigned int i = 0; i < ARRAY_SIZE(n); i++) {
        double h = (x_-_x)/(n[i]+1);            // Step length
        double* a = new double[n[i]];           // a_{i,i pm 1} elements in tridiagonal matrix
        double* b = new double[n[i]];           // a_{i,i} elements in tridiagonal matrix
        double* f = new double[n[i]];           // Sorce term

        for(unsigned int j = 0; j < n[i]; j++) {        // Initialize values
            a[j] = -1;
            b[j] = 2;
            f[j] = h*100.0*exp(-10.0*(j+1)*h);
        }
        double* f_tmp = CopyOfArray(f, n[i]);     // Backup array f

        double* c = CopyOfArray(a, n[i]);
        b[2] = 6; a[0] = 9; a[1] = -3; c[2] = 5;
        double* solution = TridiagonalSolutionGeneral(a, b, c, f, n[i]);
        double s1 = solution[0];
        double s2 = solution[1];
        double s3 = solution[2];
        double s4 = solution[3];

        delete [] a;        // Deallocate dynamic memory
        delete [] b;
        delete [] f;
    }
    return 0;
}
template<class T> T* CopyOfArray(T* arr, unsigned int n) {
    T* arr2 = new T[n];
    for(unsigned int i = 0; i < n; i++)
        arr2[i] = arr[i];
    return arr2;
}
template<class T> T* TridiagonalSolutionGeneral(T* a, T* b, T* c, T* f, unsigned int n) {
    /* a are elements a_{i,i-1} in tridiagonal matrix, a[0] = a_{2,1}
     * b are elements a_{i,i} in tridiagonal matrix, b[0] = a_{1,1}
     * c are elements a_{i,i+1} in tridiagonal matrix, c[0] = a_{1,2}
     * f are sorce elements b_{i}, f[0] = b_{1}
     */
    T temp;
    T* start = f;
    T* end = &f[n-1];
    while(f < end) {
        temp = (*a++)/(*b++);
        *b -= temp*(*c++);          // Eq (4) in report
        *f -= temp*(*f++);          // Eq (3)
    }
    while(f > start)
        *f -= (*--c)/(*b--)*(*f--); // Eq (6)
    do
        *f /= *b++;                 // Eq (7)
    while(f++ < end);
    return start;
}
