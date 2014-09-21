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
template<MatrixType Type, class T> T* GetMatrixDiagonal(Matrix<Type, T>* matrix);
template<class T> T RelativeError(T* u, T* v, unsigned int n);
template<class T> T* Sort(T* u, unsigned int n);
template<class T> void WriteArrayToFile(ofstream* file, T* array, unsigned int n);

int main() {
    unsigned int n[] = {4};
    FLOAT omega[] = {0.01, 0.5, 1, 5};          // Oscillator frequency
    FLOAT _r = 0, r_ = 1E7;                     // Solution interval
    ofstream timefile, timefile1, errorfile, errorfile1;
    timefile.open("time.dat");
    timefile1.open("time1.dat");
    errorfile.open("error.dat");
    errorfile1.open("error1.dat");

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
        errorfile << n[i] << " & ";
        timefile1 << n[i] << " & ";
        errorfile1 << n[i] << " & ";
        auto eigv = Matrix<MatrixType::Square, FLOAT>(n[i]);    // Eigenvectors
        eigv.Diagonal(0,1);

        auto matrix = Matrix<MatrixType::tqli, FLOAT>(n[i]);    // Benchmark, tqli in lib.cpp
        matrix.Clear();
        matrix.Diagonal(0, d);
        matrix.Diagonal(1,-1);
        clock_t t0 = clock();
        auto l0 = matrix.tqli(eigv);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        Sort(l0, n[i]);

        //MatrixCout(matrix);
        //MatrixCout(eigv);

        FLOAT* l1;
        auto matrix1 = Matrix<MatrixType::Square, FLOAT>(n[i]);
        if(n[i] <= 100) {
            matrix1.Clear();
            matrix1.Diagonal(0, d);
            matrix1.Diagonal(1,-1);
            matrix1.Diagonal(-1,-1);
            //MatrixCout(matrix1);
            t0 = clock();
            matrix1.JacobiMethod(1e-6);
            timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
            l1 = GetMatrixDiagonal(&matrix1);
            errorfile << RelativeError(l0, Sort(l1, n[i]), n[i]) << " & ";
            delete [] l1;
        } else {
            timefile << "- & ";
            errorfile << "- & ";
        }

        auto matrix2 = Matrix<MatrixType::Symmetric, FLOAT>(n[i]);
        if(n[i] <= 100) {
            matrix2.Clear();
            matrix2.Diagonal(0, d);
            matrix2.Diagonal(1,-1);
            //MatrixCout(matrix2);
            t0 = clock();
            auto l2 = matrix2.JacobiMethod(1e-6);
            timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
            errorfile << RelativeError(l0, Sort(l2, n[i]), n[i]) << " & ";
            //MatrixCout(matrix2);
        } else{
            timefile << "- & ";
            errorfile << "- & ";
        }

        matrix2.Clear();
        matrix2.Diagonal(0, d);
        matrix2.Diagonal(1,-1);
        t0 = clock();
        auto l3 = matrix2.JacobiMethodFD_1(1e-6);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        errorfile << RelativeError(l0, Sort(l3, n[i]), n[i]) << " & ";
        //MatrixCout(matrix2);

        matrix2.Clear();
        matrix2.Diagonal(0, d);
        matrix2.Diagonal(1,-1);
        t0 = clock();
        auto l4 = matrix2.JacobiMethodRD(1e-6);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        errorfile << RelativeError(l0, Sort(l4, n[i]), n[i]) << " & ";
        //MatrixCout(matrix2);

        matrix2.Clear();
        matrix2.Diagonal(0, d);
        matrix2.Diagonal(1,-1);
        t0 = clock();
        auto l5 = matrix2.JacobiMethodFC(1e-6);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        errorfile << RelativeError(l0, Sort(l5, n[i]), n[i]) << " & ";
        //MatrixCout(matrix2);

        matrix2.Clear();
        matrix2.Diagonal(0, d);
        matrix2.Diagonal(1,-1);
        t0 = clock();
        auto l6 = matrix2.JacobiMethodRC(1e-6);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        errorfile << RelativeError(l0, Sort(l6, n[i]), n[i]) << " & ";
        //MatrixCout(matrix6);

        matrix2.Clear();
        matrix2.Diagonal(0, d);
        matrix2.Diagonal(1,-1);
        t0 = clock();
        auto l7 = matrix2.JacobiMethodFR(1e-6);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        errorfile << RelativeError(l0, Sort(l7, n[i]), n[i]) << " & ";
        //MatrixCout(matrix7);

        matrix2.Clear();
        matrix2.Diagonal(0, d);
        matrix2.Diagonal(1,-1);
        t0 = clock();
        auto l8 = matrix2.JacobiMethodRR(1e-6);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
        errorfile << RelativeError(l0, Sort(l8, n[i]), n[i]) << " \\\\" << endl;
        //MatrixCout(matrix8);

        matrix2.Clear();
        matrix2.Diagonal(0, d);
        matrix2.Diagonal(1,-1);
        t0 = clock();
        auto l9 = matrix2.JacobiMethodFD(-4);
        timefile1 << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
        errorfile1 << RelativeError(l0, Sort(l3, n[i]), n[i]) << " \\\\" << endl;

        delete [] rho;
        delete [] d;
    }
    timefile.close();
    timefile1.close();
    errorfile.close();
    errorfile1.close();
    return 0;
}
template<class T> T* CopyOfArray(T* arr, unsigned int n) {
    T* arr2 = new T[n];
    for(unsigned int i = 0; i < n; i++)
        arr2[i] = arr[i];
    return arr2;
}
template<MatrixType Type, class T> T* GetMatrixDiagonal(Matrix<Type, T>* matrix) {
    T* d = new T[matrix->n];
    for(unsigned int i = 0; i < matrix->n; i++)
        d[i] = (*matrix)(i,i);
    return d;
}

template<class T> T* Sort(T* u, unsigned int n) {
    T temp;
    for(unsigned int i = 0; i < n-1; i++) {
        for(unsigned int j = i+1; j< n; j++) {
            if(u[i] > u[j]) {
                temp = u[i];
                u[i] = u[j];
                u[j] = temp;
            }
        }
    }
    return u;
}
template<class T> T RelativeError(T* u, T* v, unsigned int n) {
    T e, error = 0;
    for(unsigned int i = 0; i < n; i++) {
        e = abs((v[i]-u[i])/u[i]);
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
