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
template<class T> void CoutArray(T* array, unsigned int n);
template<MatrixType Type, class T> T* GetMatrixDiagonal(Matrix<Type, T>* matrix);
template<class T> T* NormalizeEigenVectors(T** v, unsigned int n);
template<class T> T RelativeError(T* u, T* v, unsigned int n);
template<class T> T* SortEigenValues(T* u, T** v, unsigned int n);
template<class T> void WriteArrayToFile(ofstream* file, T* array, unsigned int n);

int main() {
    unsigned int n[] = {4};
    FLOAT omega[] = {0.01, 0.5, 1, 5};          // Oscillator frequency
    FLOAT k = M*omega[2];
    FLOAT alpha = pow(HBAR*HBAR/(M*k),0.25);
    FLOAT _rho = 0;
    FLOAT rho_[] = {1,2,4,10,100,1000};
    ofstream timefile, errorfile, numfile;
    timefile.open("time.dat");
    errorfile.open("error.dat");
    numfile.open("numberoftransformations.dat");

    /* Compare different eigenvalue solvers */
    for(unsigned int i = 0; i < ARRAY_SIZE(n); i++) {
        unsigned int num;
        FLOAT h = (rho_[0]-_rho)/(n[i]+1);         // Step length
        FLOAT* rho = new FLOAT[n[i]+2];         // the variable rho
        FLOAT* d = new FLOAT[n[i]];             // Diagonal elements

        rho[n[i]+1] = rho_[0];
        for(unsigned int j = 0; j < n[i]; j++) {
            rho[j+1] = (j+1)*h;
            d[j] = 2 + h*h*rho[j+1]*rho[j+1];
        }
        timefile << n[i] << " & ";
        errorfile << n[i] << " & ";
        numfile << n[i] << " & ";
        auto eigv = Matrix<MatrixType::Square, FLOAT>(n[i]);    // Eigenvectors
        eigv.Clear();
        eigv.Diagonal(0,1);

        auto matrix = Matrix<MatrixType::tqli, FLOAT>(n[i]);    // Benchmark, tqli in lib.cpp
        matrix.Clear();
        matrix.Diagonal(0, d);
        matrix.Diagonal(1,-1);
        clock_t t0 = clock();
        auto l0 = matrix.tqli(eigv);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        SortEigenValues(l0, (FLOAT**)eigv, n[i]);

        auto matrix0 = Matrix<MatrixType::TridiagonalSymmetric, FLOAT>(n[i]);
        matrix0.Diagonal(0, d);
        matrix0.Diagonal(1,-1);
        t0 = clock();
        auto l = matrix0.QRalgorithm(1e-6);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        errorfile << RelativeError(l0, SortEigenValues(l, (FLOAT**)eigv, n[i]), n[i]) << " & ";

        if(n[i] <= 100) {
            auto matrix1 = Matrix<MatrixType::Square, FLOAT>(n[i]);
            matrix1.Clear();
            matrix1.Diagonal(0, d);
            matrix1.Diagonal(1,-1);
            matrix1.Diagonal(-1,-1);
            t0 = clock();
            matrix1.JacobiMethod(1e-6, num);
            timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
            l = GetMatrixDiagonal(&matrix1);
            errorfile << RelativeError(l0, SortEigenValues(l, (FLOAT**)eigv, n[i]), n[i]) << " & ";
            numfile << num << " & ";
            delete [] l;
        } else {
            timefile << "- & ";
            errorfile << "- & ";
            numfile << "- & ";
        }

        auto matrix2 = Matrix<MatrixType::Symmetric, FLOAT>(n[i]);
        if(n[i] <= 100) {
            matrix2.Clear();
            matrix2.Diagonal(0, d);
            matrix2.Diagonal(1,-1);
            t0 = clock();
            l = matrix2.JacobiMethod(1e-6, num);
            timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
            errorfile << RelativeError(l0, SortEigenValues(l, (FLOAT**)eigv, n[i]), n[i]) << " & ";
            numfile << num << " & ";
        } else {
            timefile << "- & ";
            errorfile << "- & ";
            numfile << "- & ";
        }

        matrix2.Clear();
        matrix2.Diagonal(0, d);
        matrix2.Diagonal(1,-1);
        t0 = clock();
        l = matrix2.JacobiMethodFD(1e-6, num);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        errorfile << RelativeError(l0, SortEigenValues(l, (FLOAT**)eigv, n[i]), n[i]) << " & ";
        numfile << num << " & ";

        matrix2.Clear();
        matrix2.Diagonal(0, d);
        matrix2.Diagonal(1,-1);
        t0 = clock();
        l = matrix2.JacobiMethodRD(1e-6, num);
        timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
        errorfile << RelativeError(l0, SortEigenValues(l, (FLOAT**)eigv, n[i]), n[i]) << " & ";
        numfile << num << " & ";

        if(n[i] <= 100) {
            matrix2.Clear();
            matrix2.Diagonal(0, d);
            matrix2.Diagonal(1,-1);
            t0 = clock();
            l = matrix2.JacobiMethodFC(1e-6, num);
            timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
            errorfile << RelativeError(l0, SortEigenValues(l, (FLOAT**)eigv, n[i]), n[i]) << " & ";
            numfile << num << " & ";
        } else {
            timefile << "- & ";
            errorfile << "- & ";
            numfile << "- & ";
        }

        if(n[i] <= 100) {
            matrix2.Clear();
            matrix2.Diagonal(0, d);
            matrix2.Diagonal(1,-1);
            t0 = clock();
            l = matrix2.JacobiMethodRC(1e-6, num);
            timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
            errorfile << RelativeError(l0, SortEigenValues(l, (FLOAT**)eigv, n[i]), n[i]) << " & ";
            numfile << num << " & ";
        } else {
            timefile << "- & ";
            errorfile << "- & ";
            numfile << "- & ";
        }

        if(n[i] <= 100) {
            matrix2.Clear();
            matrix2.Diagonal(0, d);
            matrix2.Diagonal(1,-1);
            t0 = clock();
            l = matrix2.JacobiMethodFR(1e-6, num);
            timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " & ";
            errorfile << RelativeError(l0, SortEigenValues(l, (FLOAT**)eigv, n[i]), n[i]) << " & ";
            numfile << num << " & ";
        } else {
            timefile << "- & ";
            errorfile << "- & ";
            numfile << "- & ";
        }

        if(n[i] <= 100) {
            matrix2.Clear();
            matrix2.Diagonal(0, d);
            matrix2.Diagonal(1,-1);
            t0 = clock();
            l = matrix2.JacobiMethodRR(1e-6, num);
            timefile << (double)(clock()-t0)/CLOCKS_PER_SEC << " \\\\" << endl;
            errorfile << RelativeError(l0, SortEigenValues(l, (FLOAT**)eigv, n[i]), n[i]) << " \\\\" << endl;
            numfile << num << " \\\\" << endl;
        } else {
            timefile << "- \\\\" << endl;
            errorfile << "- \\\\" << endl;
            numfile << "- \\\\" << endl;
        }

        delete [] rho;
        delete [] d;
    }
    timefile.close();
    errorfile.close();
    numfile.close();

    /* Single electron harmonic oscillator */
    unsigned int nstep[] = {/*3,*/4,10,20,50,100,1000};
    FLOAT L1[] = {3,7,11};
    errorfile.open("error2header.dat");
    errorfile << "$n_{\\text{step}} \\textbackslash \\rho_{\\text{max}}$";
    for(unsigned int i = 0; i < ARRAY_SIZE(rho_); i++)
        errorfile << " & " << rho_[i];
    errorfile << " \\\\";
    errorfile.close();
    errorfile.open("error2.dat");
    for(unsigned int i = 0; i < ARRAY_SIZE(nstep); i++) {
        errorfile << nstep[i];
        for(unsigned int j = 0; j < ARRAY_SIZE(rho_); j++) {
            FLOAT h = (rho_[j]-_rho)/(nstep[i]+1);         // Step length
            FLOAT* rho = new FLOAT[nstep[i]+2];         // the variable rho
            FLOAT* d = new FLOAT[nstep[i]];             // Diagonal elements

            rho[nstep[i]+1] = rho_[j];
            for(unsigned int k = 0; k < nstep[i]; k++) {
                rho[k+1] = (k+1)*h;
                d[k] = 2+h*h*rho[k+1]*rho[k+1];
            }

            auto eigv = Matrix<MatrixType::Square, FLOAT>(nstep[i]);    // Eigenvectors
            eigv.Clear();
            eigv.Diagonal(0,1);

            auto matrix = Matrix<MatrixType::tqli, FLOAT>(nstep[i]);    // Benchmark, tqli in lib.cpp
            matrix.Clear();
            matrix.Diagonal(0, d);
            matrix.Diagonal(1,-1);
            auto l = matrix.tqli(eigv);
            SortEigenValues(l, (FLOAT**)eigv, nstep[i]);
            for(unsigned int k = 0; k < ARRAY_SIZE(L1); k++)
                l[k] /= h*h;
            errorfile << " & " << RelativeError(L1, l, ARRAY_SIZE(L1));

            delete [] rho;
            delete[] d;
        }
        errorfile << " \\\\" << endl;
    }
    errorfile.close();

    return 0;
}
template<class T> T* CopyOfArray(T* arr, unsigned int n) {
    T* arr2 = new T[n];
    for(unsigned int i = 0; i < n; i++)
        arr2[i] = arr[i];
    return arr2;
}
template<class T> void CoutArray(T* array, unsigned int n) {
    for(unsigned int i = 0; i < n; i++)
        cout << array[i] << '\t';
    cout << endl;
}
template<MatrixType Type, class T> T* GetMatrixDiagonal(Matrix<Type, T>* matrix) {
    T* d = new T[matrix->n];
    for(unsigned int i = 0; i < matrix->n; i++)
        d[i] = (*matrix)(i,i);
    return d;
}
template<class T> T* NormalizeEigenVectors(T** v, unsigned int n) {
    T* f = new T[n];
    for(unsigned int i = 0; i < n; i++) {
        f[i] = 0;
        for(unsigned int j = 0; j < n; j++)
            f[i] += v[j][i]*v[j][i];
        f[i] = sqrt(f[i]);
    }
}
template<class T> T* SortEigenValues(T* u, T** v, unsigned int n) {
    T temp;
    for(unsigned int i = 0; i < n-1; i++) {
        for(unsigned int j = i+1; j< n; j++) {
            if(u[i] > u[j]) {
                temp = u[i];
                u[i] = u[j];
                u[j] = temp;
                for(unsigned int k = 0; k < n; k++) {
                    temp = v[k][i];
                    v[k][i] = v[k][j];
                    v[k][j] = temp;
                }
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
