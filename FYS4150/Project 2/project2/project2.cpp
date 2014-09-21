/*
 *  FYS4150 - Computational Physics - Project 2
 *
 *  Written by: Eimund Smestad
 *
 *  16.09.2014
 *
 *  c11 compiler
 */

// Include files
#include <cmath>
#include <iostream>
#include <fstream>
#include "time.h"

// My include files
#include "Matrix.h"

// Namespaces
using namespace std;

// Constants
#define HBAR    6.62606957E-34                  // Plancks constant
#define M       9.10938291E-31                  // Electron mass

// Macros
#define FLOAT double                            // Choosen floating point precision
#define ARRAY_SIZE(X) sizeof(X)/sizeof(*X)      // Calculate static array size at compile time

// Declartion of global functions
template<class T> T* CopyOfArray(T* arr, unsigned int n);
//template<MatrixType Type, class T> void CoutMatrix(Matrix<Type, T>& matrix);
template<class T> T RelativeError(T* u, T* v, unsigned int n);
template<class T> void WriteArrayToFile(ofstream* file, T* array, unsigned int n);

int main() {
    unsigned int n[] = {1000};
    FLOAT omega[] = {0.01, 0.5, 1, 5};          // Oscillator frequency
    FLOAT _r = 0, r_ = 1E7;                     // Solution interval
    ofstream timefile;
    timefile.open("time.dat");

    for(unsigned int i = 0; i < ARRAY_SIZE(n); i++) {

        FLOAT k = M*omega[2];
        FLOAT alpha = pow(HBAR*HBAR/(M*k),0.25);
        FLOAT _rho = _r/alpha;
        FLOAT rho_ = r_/alpha;
        FLOAT h = (rho_-_rho)/(n[i]+1);         // Step length
        FLOAT* rho = new FLOAT[n[i]+2];         // the variable rho
        FLOAT* d = new FLOAT[n[i]];             // Diagonal elements

        rho[n[i]+1] = rho_;
        for(unsigned int j = 0; j < n[i]; j++) {
            rho[j+1] = (j+1)*h;
            d[j] = 2 + 0.5*k*pow(h*alpha*rho[j+1],2);
        }
        timefile << n[i] << " & ";
        auto eigv = Matrix<MatrixType::Square, FLOAT>(n[i]);    // Eigenvectors
        eigv.Diagonal(0,1);

        auto matrix = Matrix<MatrixType::tqli, FLOAT>(n[i]);    // Benchmark, tqli in lib.cpp
        matrix.Diagonal(0, d);
        matrix.Diagonal(1,-1);
        clock_t t0 = clock();
        matrix.tqli(eigv);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        //MatrixCout(matrix);
        //MatrixCout(eigv);

        if(n[i] <= 100) {
            auto matrix1 = Matrix<MatrixType::Square, FLOAT>(n[i]);
            matrix1.Diagonal(0, d);
            matrix1.Diagonal(1,-1);
            matrix1.Diagonal(-1,-1);
            //MatrixCout(matrix1);
            t0 = clock();
            matrix1.JacobiMethod(1e-3,1);
            timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        } else
            timefile << "- & ";

        if(n[i] <= 100) {
            auto matrix2 = Matrix<MatrixType::Symmetric, FLOAT>(n[i]);
            matrix2.Diagonal(0, d);
            matrix2.Diagonal(1,-1);
            //MatrixCout(matrix2);
            t0 = clock();
            matrix2.JacobiMethod(1e-3,1);
            timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        } else
            timefile << "- & ";

        auto matrix3 = Matrix<MatrixType::Symmetric, FLOAT>(n[i]);
        matrix3.Diagonal(0, d);
        matrix3.Diagonal(1,-1);
        t0 = clock();
        matrix3.JacobiMethodFaster(1e-3,1);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & " << endl;

        auto matrix4 = Matrix<MatrixType::Symmetric, FLOAT>(n[i]);
        matrix4.Diagonal(0, d);
        matrix4.Diagonal(1,-1);
        t0 = clock();
        matrix4.JacobiMethodFastest(1e-3,1);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;

        auto matrix5 = Matrix<MatrixType::Symmetric, FLOAT>(5);
        matrix5(0,0) = 1;
        matrix5(0,1) = 2;
        matrix5(0,2) = 3;
        matrix5(0,3) = 4;
        matrix5(0,4) = 5;
        matrix5(1,1) = 6;
        matrix5(1,2) = 7;
        matrix5(1,3) = 8;
        matrix5(1,4) = 9;
        matrix5(2,2) = 10;
        matrix5(2,3) = 11;
        matrix5(2,4) = 12;
        matrix5(3,3) = 13;
        matrix5(3,4) = 14;
        matrix5(4,4) = 15;
        //MatrixCout(matrix4);

        delete [] rho;
        delete [] d;
    }
    timefile.close();
    return 0;
}
template<class T> T* CopyOfArray(T* arr, unsigned int n) {
    T* arr2 = new T[n];
    for(unsigned int i = 0; i < n; i++)
        arr2[i] = arr[i];
    return arr2;
}
template<class T> T RelativeError(T* u, T* v, unsigned int n) {
    T e, error = 0;
    for(unsigned int i = 0; i < n; i++) {
        e = abs((v[0]-u[0])/u[0]);
        if(e > error)
            error = e;
    }
    return log10(error);
}
template<class T> void WriteArrayToFile(ofstream* file, T* array, unsigned int n) {
    for(unsigned int i = 0; i < n-1; i++)
        *file << array[i] << '\t';
    *file << array[n-1];
}
