/*
 *  FYS4150 - Computational Physics - Project 1
 *
 *  Written by: Eimund Smestad
 *
 *  30.08.2014
 */

// Include files
#include <cmath>
#include <iostream>
#include "armadillo"

// Namespaces
using namespace arma;
using namespace std;

// Macros
#define ARRAY_SIZE(X) sizeof(X)/sizeof(*X)      // Calculate static array size at compile time

// Classes
template<class T> class GetPointer {
};
template<class T> class GetPointer<T*> {
    typedef T Type;                             // Used to get the type the pointer points to at compile time
};
template<unsigned int Type, class V> class TridiagonalMatrix {              // General tridiagonal matrix solver
    public: static void Solve(V* a, V* b, V* c, V* f, unsigned int n) {     // 8n FLOPS
        /* a are elements a_{i,i-1} in tridiagonal matrix, a[0] = a_{2,1}
         * b are elements a_{i,i} in tridiagonal matrix, b[0] = a_{1,1}
         * c are elements a_{i,i+1} in tridiagonal matrix, c[0] = a_{1,2}
         * f are sorce elements b_{i}, f[0] = b_{1}
         */

        V temp;
        V* start = f;
        V* end = &f[n-1];
        while(f < end) {
            temp = (*a++)/(*b++);
            *b -= temp*(*c++);          // Eq (4) in report
            *f -= temp*(*f++);          // Eq (3)
        }
        while(f > start) {
            *f /= (*b--);               // Eq (7)
            *f -= (*--c)*(*f--);        // Eq (6)
        }
        *f /= *b;

        /* Testet matrise
         *
         *  2  -1   0   0   2.70671              0.0392731
         *  9   2  -1   0   0.366313     =>     -1.92124
         *  0  -3   6   5   0.049575            -0.674227
         *  0   0  -1   2   0.00670925          -0.333759
         *
         *  2  -1   0   0   2.70671              2.40633
         * -1  2   -1   0   0.366313     =>      2.10595
         *  0  -1   2  -1   0.049575             1.43925
         *  0   0  -1   2   0.00670925           0.72298
         */
    }
};
template<class V> class TridiagonalMatrix<1,V> {                            // Tridiagonal matrix -1,2,-1 without memory usage
    public: static void Solve(V* f, unsigned int n) {                       // 6n FLOPS
        Solve(f, n, n);
    }
    public: static void Solve(V* f, unsigned int n, unsigned int cutoff) {
        /*
         * f are sorce elements b_{i}, f[0] = b_{1}
         * The values of i >= cutoff is when (i+1)/i is approximated to 1
         */

        V* start = f;
        V* mid = &f[cutoff-1];
        V* end = &f[n-1];
        unsigned int i = 1;
        while(f < mid)
            *f += (V)i++/i*(*f++);      // Eq (8)
        while(f < end)
            *f += *f++;                 // Eq (8) with (i+1)/i = 1
        i++;
        while(f > mid)
            *f += *f--;                 // Eq (10) with (i+1)/i = 1
        while(f > start) {
            *f /= (V)i--/i;             // Eq (11)
            *f += *f--;                 // Eq (10)
        }
        *f /= 2;
    }
};

// Declartion of global functions
template<class T> T* CopyOfArray(T* arr, unsigned int n);

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

        TridiagonalMatrix<0,GetPointer<decltype(f)>::Type>::Solve(a, b, a, f, n[i]);
        double s1 = f[0];
        double s2 = f[1];
        double s3 = f[2];
        double s4 = f[3];

        delete [] f;
        f = f_tmp;
        f_tmp = CopyOfArray(f, n[i]);
        TridiagonalMatrix<1,GetPointer<decltype(f)>::Type>::Solve(f, n[i],3);
        double t1 = f[0];
        double t2 = f[1];
        double t3 = f[2];
        double t4 = f[3];

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
