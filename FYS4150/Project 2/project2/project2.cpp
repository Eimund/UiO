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
template<class T> T RelativeError(T* u, T* v, unsigned int n);
template<class T> void WriteArrayToFile(ofstream* file, T* array, unsigned int n);

int main() {
    unsigned int n[] = {4};
    FLOAT omega[] = {0.01, 0.5, 1, 5};          // Oscillator frequency
    FLOAT _r = 0, r_ = 1E7;                       // Solution interval

    for(unsigned int i = 0; i < ARRAY_SIZE(n); i++) {

        FLOAT k = M*omega[2];
        FLOAT alpha = pow(HBAR*HBAR/(M*k),0.25);
        FLOAT _rho = _r/alpha;
        FLOAT rho_ = r_/alpha;
        FLOAT h = (rho_-_rho)/(n[i]+1);         // Step length
        FLOAT* rho = new FLOAT[n[i]+2];         // the variable rho
        FLOAT* d = new FLOAT[n[i]];             // Diagonal elements
        FLOAT* e = new FLOAT[n[i]];             // Off diagonal elements
        FLOAT** z = new FLOAT*[n[i]];           // Eigenvector;

        rho[n[i]+1] = rho_;
        for(unsigned int j = 0; j < n[i]; j++) {
            rho[j+1] = (j+1)*h;
            d[j] = 2 + 0.5*k*pow(h*alpha*rho[j+1],2);
            e[j] = -1;
            z[j] = new FLOAT[n[i]];
        }
        auto matrix = Matrix<MatrixType::Symmetric, FLOAT>(n[i]);
        matrix.Diagonal(0, d);
        matrix.Diagonal(1, -1);

        FLOAT p1 = d[0];
        FLOAT p2 = d[1];
        FLOAT p3 = d[2];
        FLOAT p4 = d[3];
        tqli(d, e, n[i], z);
        FLOAT d1 = d[0];
        FLOAT d2 = d[1];
        FLOAT d3 = d[2];
        FLOAT d4 = d[3];

        delete [] rho;
        delete [] d;
        delete [] e;
        for(unsigned int j = 0; j < n[i]; j++)
            delete [] z[j];
        delete [] z;
    }
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
