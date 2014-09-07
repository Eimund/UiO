/*
 *  FYS4150 - Computational Physics - Project 1
 *
 *  Written by: Eimund Smestad
 *
 *  30.08.2014
 *
 *  c11 compiler
 */

// Include files
#include <cmath>
#include <fstream>
#include "armadillo"
#include "time.h"

// Namespaces
using namespace arma;
using namespace std;

// Macros
#define FLOAT double
#define ARRAY_SIZE(X) sizeof(X)/sizeof(*X)      // Calculate static array size at compile time

// Classes
template<class T> class GetPointer {
};
template<class T> class GetPointer<T*> {
    public: typedef T Type;                             // Used to get the type the pointer points to at compile time
};
template<unsigned int Type, class T> class TridiagonalMatrix {              // General tridiagonal matrix solver
    public: static T* Solve(T* a, T* b, T* c, T* f, unsigned int n) {     // 8n FLOPS
        /* a are elements a_{i,i-1} in tridiagonal matrix, a[0] = a_{2,1}
         * b are elements a_{i,i} in tridiagonal matrix, b[0] = a_{1,1}
         * c are elements a_{i,i+1} in tridiagonal matrix, c[0] = a_{1,2}
         * f are source elements b_{i}, f[0] = b_{1}
         */

        T temp;
        T* start = f;
        T* end = &f[n-1];
        while(f != end) {
            temp = (*a++)/(*b++);
            *b -= temp*(*c++);          // Eq (4) in report
            *f -= temp*(*f++);          // Eq (3)
        }
        while(f != start) {
            *f /= (*b--);               // Eq (7)
            *f -= (*--c)*(*f--);        // Eq (6)
        }
        *f /= *b;
        return f;

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
template<class T> class TridiagonalMatrix<1,T> {                            // Tridiagonal matrix -1,2,-1 without memory usage
    public: static T* Solve(T* f, unsigned int n) {                         // 6n FLOPS => 2n FLOPS when n >> cutoff
        return Solve(f, n, n);
    }
    public: static T* Solve(T* f, unsigned int n, unsigned int cutoff) {
        /*
         * f are source elements b_{i}, f[0] = b_{1}
         * The values of i >= cutoff is when (i+1)/i is approximated to 1
         */

        T* start = f;
        T* mid = &f[n < cutoff ? n-1 : cutoff-1];
        T* end = &f[n-1];
        T i = 1;                        // Faster to use double than int here
        while(f != mid)
            *f += i++/i*(*f++);         // Eq (8)
        while(f != end)
            *f += *f++;                 // Eq (8) with (i+1)/i = 1
        i++;
        while(f != mid)
            *f += *f--;                 // Eq (10) with (i+1)/i = 1
        while(f != start) {
            *f *= (T)1/(i--)*i;         // Eq (11) (faster than *f /= i--/i)
            *f += *f--;                 // Eq (10)
        }
        *f /= (T)2;
        return f;
    }
};
template<class T> class TridiagonalMatrix<2,T> {                            // Tridiagonal matrix -1,2,-1 with memory usage
    public: T* factor;                   // Precalculated values
    private: unsigned int _cutoff;
    public: class _cutoff_ {             // property class
        private: TridiagonalMatrix<2,T>* owner;
        public: _cutoff_(TridiagonalMatrix<2,T>* owner) : owner(owner) {
            this->owner->_cutoff = 0;
        }
        public: unsigned int & operator = (const unsigned int& cutoff) { // set function
            if(owner->_cutoff != cutoff) {
                T* arr = new T[cutoff];
                T* start;
                if(owner->_cutoff < cutoff) {       // Extend values
                    start = &arr[owner->_cutoff]-1;
                    arr = &arr[cutoff-1];
                    double i = cutoff+1;
                    while(arr != start)
                        *arr-- = (T)1/(i--)*i;          // Calculate (i-1)/i values
                }
                if(owner->_cutoff!= 0) {
                    start = owner->factor-1;
                    if(owner->_cutoff < cutoff)
                        owner->factor= &owner->factor[owner->_cutoff-1];
                    else {
                        arr = &arr[cutoff-1];
                        owner->factor= &owner->factor[cutoff-1];
                    }
                    while(owner->factor != start)   // Copy already calculated values
                        *arr-- = *(owner->factor--);
                    delete [] ++start;
                }
                owner->factor = ++arr;
            }
            return owner->_cutoff= cutoff;
        }
        public: operator unsigned int () const {    // get function
            return owner->_cutoff;
        }
    } cutoff;   // Number of precalculated values
    public: TridiagonalMatrix(unsigned int cutoff) : cutoff(this) {
        this->cutoff = cutoff;
    }
    public: ~TridiagonalMatrix() {
        delete [] factor;
    }
    public: T* Solve(T* f, unsigned int n) {         // 4n FLOPS => 2n FLOPS for n >> cutoff
        /*
         * f are source elements b_{i}, f[0] = b_{1}
         */

        T* start = f;
        T* mid = &f[n < _cutoff ? n-1 : _cutoff-1];
        T* end = &f[n-1];
        T* fac = factor;
        while(f != mid)
            *f += (*f++)*(*fac++);      // Eq (8)
        while(f != end)
            *f += *f++;                 // Eq (8) with (i+1)/i = 1
        while(f != mid)
            *f += *f--;                 // Eq (10) with (i+1)/i = 1
        while(f != start) {
            *f *= *fac--;               // Eq (11)
            *f += *f--;                 // Eq (10)
        }
        *f /= (T)2;
        return f;
    }
};
template<class T> class LUMatrix {
    public: inline static void Solve(Mat<T>* X, T* f) {
        Mat<T> L, U;
        lu(L, U, *X);                   // Armadillo LU decomposition
        int i = 0, j, n = X->n_rows;
        T* f_;
        while(++i < n) {                // Forward solve Ly = f
            f_ = &f[i];
            for(j = 0; j < i; j++)
                *f_ -= L(i,j)*f[j];
        }
        n--;
        while(i--) {                    // Backward solve Ux = y
            f_ = &f[i];
            for(j = n; j > i; --j)
                *f_ -= U(i,j)*f[j];
            *f_ /= U(i,i);
        }
    }
};

// Declartion of global functions
template<class T> T* CopyOfArray(T* arr, unsigned int n);
template<class T> void WriteArrayToFile(ofstream* file, T* array, unsigned int n);

int main() {
    unsigned int n[] = {4,100,1000,10000,100000000};  // Size of matrix to solve
    double _x = 0, x_ = 1;                      // Solution interval
    char filename[50];
    ofstream timefile, file;
    timefile.open("time.dat");

    for(unsigned int i = 0; i < ARRAY_SIZE(n); i++) {
        mat A;                                  // Matrix for armadillo
        double h = (x_-_x)/(n[i]+1);            // Step length
        double* a = new double[n[i]];           // a_{i,i pm 1} elements in tridiagonal matrix
        double* b = new double[n[i]];           // a_{i,i} elements in tridiagonal matrix
        double* f = new double[n[i]];           // Sorce term
        double* x = new double[n[i]];           // the variable x
        double* u = new double[n[i]];           // the closed form solution

        if(n[i] <= 1000) {
            A = mat(n[i],n[i],fill::zeros);

            sprintf(filename,"result_\%d.dat",n[i]);
            file.open(filename);
            WriteArrayToFile(&file, x, n[i]);
            file << endl;
            WriteArrayToFile(&file, u, n[i]);
            file << endl;
        }
        for(unsigned int j = 0; j < n[i]; j++) {        // Initialize values
            a[j] = -1;
            b[j] = 2;
            x[j] = (j+1)*h;
            f[j] = h*h*100.0*exp(-10.0*x[j]);
            u[j] = (1-(1-exp(-10))*x[j]-exp(-10*x[j]));

            if(n[i] <= 1000) {
                A(j,j)= 2;
                if(j) {
                    A(j,j-1) = -1;
                    A(j-1,j) = -1;
                }
            }
        }        
        double* f_tmp = CopyOfArray(f, n[i]);     // Backup array f
        timefile << n[i] << " & ";

        clock_t t0 = clock();
        f = TridiagonalMatrix<0,GetPointer<decltype(f)>::Type>::Solve(a, b, a, f, n[i]); // Solve general tridiagonal matrix
        auto t1 = (double)(clock()-t0)/CLOCKS_PER_SEC;
        timefile << t1 << " & ";
        if(n[i] <= 1000) {
            WriteArrayToFile(&file, f, n[i]);
            file << endl;
        }

        delete [] f;
        f = f_tmp;
        f_tmp = CopyOfArray(f, n[i]);
        t0 = clock();
        f = TridiagonalMatrix<1,GetPointer<decltype(f)>::Type>::Solve(f, n[i]);         // Solve -1,2,-1 tridiagonal matrix without memory usage
        auto t2 = (double)(clock()-t0)/CLOCKS_PER_SEC;
        timefile << t2 << " & ";
        if(n[i] <= 1000) {
            WriteArrayToFile(&file, f, n[i]);
            file << endl;
        }

        delete [] f;
        f = f_tmp;
        f_tmp = CopyOfArray(f, n[i]);
        t0 = clock();
        f = TridiagonalMatrix<2,GetPointer<decltype(f)>::Type>(n[i]).Solve(f, n[i]);    // Solve -1,2,-1 tridiagonal matrix without precalculated values but with memory usage
        auto t3 = (double)(clock()-t0)/CLOCKS_PER_SEC;
        timefile << t3 << " & ";
        if(n[i] <= 1000) {
            WriteArrayToFile(&file, f, n[i]);
            file << endl;
        }

        delete [] f;
        f = f_tmp;
        f_tmp = CopyOfArray(f, n[i]);
        auto tri = TridiagonalMatrix<2,GetPointer<decltype(f)>::Type>(n[i]);            // Precalculate values for
        t0 = clock();
        f = tri.Solve(f, n[i]);                                                         // solving -1,2,-1 tridiagonal matrix with memory usage
        auto t4 = (double)(clock()-t0)/CLOCKS_PER_SEC;
        timefile << t4 << " & ";
        if(n[i] <= 1000)
            WriteArrayToFile(&file, f, n[i]);
        double r1 = f[0];
        double r2 = f[1];
        double r3 = f[2];
        double r4 = f[3];

        delete [] f;
        f = f_tmp;
        if(n[i] <= 1000) {
            t0 = clock();
            LUMatrix<double>::Solve(&A, f);                                             // Solve with LU decomposition amardillo
            auto t5 = (double)(clock()-t0)/CLOCKS_PER_SEC;
            timefile << t5 << " \\\\" << endl;
        } else
            timefile << "- \\\\" << endl;

        double s1 = f[0];
        double s2 = f[1];
        double s3 = f[2];
        double s4 = f[3];

        delete [] a;        // Deallocate dynamic memory
        delete [] b;
        delete [] f;
        delete [] x;
        delete [] u;
        file.close();
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
template<class T> void WriteArrayToFile(ofstream* file, T* array, unsigned int n) {
    for(unsigned int i = 0; i < n-1; i++)
        *file << array[i] << '\t';
    *file << array[n-1];
}
