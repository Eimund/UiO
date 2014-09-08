/*
 *  FYS4150 - Computational Physics - Project 1
 *
 *  Written by: Eimund Smestad
 *
 *  07.09.2014
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

// Macros
#define FLOAT double                            // Choosen floating point precision
#define ARRAY_SIZE(X) sizeof(X)/sizeof(*X)      // Calculate static array size at compile time

// Declartion of global functions
template<class T> T* CopyOfArray(T* arr, unsigned int n);
template<class T> T RelativeError(T* u, T* v, unsigned int n);
template<class T> void WriteArrayToFile(ofstream* file, T* array, unsigned int n);

int main() {
    unsigned int n[] = {10,100,1000,10000,100000000};  // Size of matrix to solve
    FLOAT _x = 0, x_ = 1;                      // Solution interval
    char filename[50];
    ofstream timefile, file, errorfile;
    timefile.open("time.dat");
    errorfile.open("error.dat");

    for(unsigned int i = 0; i < ARRAY_SIZE(n); i++) {
        FLOAT h = (x_-_x)/(n[i]+1);             // Step length
        FLOAT* f = new FLOAT[n[i]];             // Sorce term
        FLOAT* x = new FLOAT[n[i]+2];           // the variable x
        FLOAT* u = new FLOAT[n[i]+2];           // the closed form solution

        x[n[i]+1] = x_;
        for(unsigned int j = 0; j < n[i]; j++) {    // Initialize values
            x[j+1] = (j+1)*h;
            f[j] = h*h*100.0*exp(-10.0*x[j+1]);
            u[j+1] = (1-(1-exp(-10))*x[j+1]-exp(-10*x[j+1]));
        }
        if(n[i] <= 1000) {
            sprintf(filename,"result_\%d.dat",n[i]);
            file.open(filename);
            WriteArrayToFile(&file, x, n[i]+2);
            file << endl;
            WriteArrayToFile(&file, u, n[i]+2);
            file << endl;
        }
        FLOAT* f_tmp = CopyOfArray(f, n[i]);       // Backup array f
        timefile << n[i] << " & ";
        errorfile << n[i] << " & ";

        auto matrix = Matrix<MatrixType::TridiagonalGeneral, FLOAT>(n[i]);
        matrix.Diagonal(0,2);
        matrix.Diagonal(-1,-1);
        matrix.Diagonal(1,-1);
        clock_t t0 = clock();
        matrix.Solve(f, n[i]);                              // Solve the general tridiagonal matrix
        auto t1 = (double)(clock()-t0)/CLOCKS_PER_SEC;
        timefile << t1 << " & ";
        errorfile << RelativeError(&u[1], f, n[i]) << " & ";

        delete [] f;
        f = f_tmp;
        f_tmp = CopyOfArray(f, n[i]);
        t0 = clock();
        Matrix<MatrixType::Tridiagonal_minus1_2_minus1_6n, FLOAT>::Solve(f,n[i]);       // Solve -1,2,-1 tridiagonal matrix without memory usage
        auto t2 = (double)(clock()-t0)/CLOCKS_PER_SEC;
        timefile << t2 << " & ";
        errorfile << RelativeError(&u[1], f, n[i]) << " & ";

        delete [] f;
        f = f_tmp;
        f_tmp = CopyOfArray(f, n[i]);
        t0 = clock();
        Matrix<MatrixType::Tridiagonal_minus1_2_minus1_6n, FLOAT>::Solve(f,n[i], 10000);// Solve -1,2,-1 tridiagonal matrix without memory usage with cutoff
        auto t3 = (double)(clock()-t0)/CLOCKS_PER_SEC;
        timefile << t3 << " & ";
        errorfile << RelativeError(&u[1], f, n[i]) << " & ";

        delete [] f;
        f = f_tmp;
        f_tmp = CopyOfArray(f, n[i]);
        t0 = clock();
        Matrix<MatrixType::Tridiagonal_minus1_2_minus1_6n, FLOAT>::SolveInt(f,n[i]);    // Solve -1,2,-1 tridiagonal matrix without memory usage, int factor
        auto t4 = (double)(clock()-t0)/CLOCKS_PER_SEC;
        timefile << t4 << " & ";
        errorfile << RelativeError(&u[1], f, n[i]) << " & ";

        delete [] f;
        f = f_tmp;
        f_tmp = CopyOfArray(f, n[i]);
        auto matrix3 = Matrix<MatrixType::Tridiagonal_minus1_2_minus1_6n, FLOAT>(n[i]);
        t0 = clock();
        matrix3.SolveTrue(f,n[i]);                                                      // Solve -1,2,-1 tridiagonal matrix without memory usage, true 6n FLOPS
        auto t5 = (double)(clock()-t0)/CLOCKS_PER_SEC;
        timefile << t5 << " & ";
        errorfile << RelativeError(&u[1], f, n[i]) << " & ";

        delete [] f;
        f = f_tmp;
        f_tmp = CopyOfArray(f, n[i]);
        t0 = clock();
        Matrix<MatrixType::Tridiagonal_minus1_2_minus1_4n, FLOAT>(n[i]).Solve(f, n[i]); // Solve -1,2,-1 tridiagonal matrix without precalculated values but with memory usage
        auto t6 = (double)(clock()-t0)/CLOCKS_PER_SEC;
        timefile << t6 << " & ";
        errorfile << RelativeError(&u[1], f, n[i]) << " & ";

        delete [] f;
        f = f_tmp;
        f_tmp = CopyOfArray(f, n[i]);
        auto tri = Matrix<MatrixType::Tridiagonal_minus1_2_minus1_4n, FLOAT>(n[i]);     // Precalculate values for
        t0 = clock();
        tri.Solve(f, n[i]);                                                             // solving -1,2,-1 tridiagonal matrix with memory usage
        auto t7 = (double)(clock()-t0)/CLOCKS_PER_SEC;
        timefile << t7 << " & ";
        errorfile << RelativeError(&u[1], f, n[i]) << " & ";
        if(n[i] <= 1000) {
            file << '0' << '\t';
            WriteArrayToFile(&file, f, n[i]);
            file << '\t' << 0;
        }

        delete [] f;
        f = f_tmp;
        if(n[i] <= 1000) {
            auto matrix2 = Matrix<MatrixType::LU_decomposition, FLOAT>(n[i]);
            matrix2.Diagonal(0,2);
            matrix2.Diagonal(-1,-1);
            matrix2.Diagonal(1,-1);
            matrix2.LU();
            t0 = clock();
            matrix2.Solve(f, n[i]);                                                     // Solve with LU decomposition amardillo
            auto t8 = (double)(clock()-t0)/CLOCKS_PER_SEC;
            timefile << t8 << " \\\\" << endl;
            errorfile << RelativeError(&u[1], f, n[i]) << " \\\\" << endl;
        } else {
            timefile << "- \\\\" << endl;
            errorfile << "- \\\\" << endl;
        }

        delete [] f;            // Deallocate dynamic memory
        delete [] x;
        delete [] u;
        file.close();
    }
    timefile.close();
    errorfile.close();
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

