/*
 *  FYS4150 - Computational Physics - Project 1
 *
 *  Written by: Eimund Smestad
 *
 *  07.09.2014
 *
 *  c11 compiler
 */

#ifndef MATRIX_H
#define MATRIX_H

#include "armadillo"

using namespace arma;

enum class MatrixType {
    TridiagonalGeneral,
    Tridiagonal_minus1_2_minus1_6n,
    Tridiagonal_minus1_2_minus1_4n,
    LU_decomposition
};

template<MatrixType Type, class T> class Matrix;
template<class M, class T> class MatrixDiagonal {
    private: M* owner;
    protected: MatrixDiagonal(M* owner) : owner(owner) {        // Fill a diagonal in a matrix with on value
    }
    public: inline void Diagonal(const int diagonal, const T value) {
        Diagonal(diagonal, value, owner->n);
    }
    public: inline void Diagonal(const int diagonal, const T value, const unsigned int n) {
        unsigned int d = abs(diagonal);
        owner->n = n;
        unsigned int nmax = n - d;
        if(diagonal > 0) {
            for(unsigned int i = 0; i < nmax; i++)
                owner->operator()(i,i+d) = value;
        } else {
            for(unsigned int i = 0; i < nmax; i++)
                owner->operator()(i+d,i) = value;
        }
    }
};
template<class T> class Matrix<MatrixType::TridiagonalGeneral, T> : MatrixDiagonal<Matrix<MatrixType::TridiagonalGeneral, T>, T> {
    private: T other;
    private: T *a, *b, *c;
    private: unsigned int _n;
    private: class _n_ {             // property class
        private: Matrix<MatrixType::TridiagonalGeneral, T>* owner;
        public: _n_(Matrix<MatrixType::TridiagonalGeneral, T>* owner) : owner(owner) {
            this->owner->_n = 1;
            owner->a = new T[0];
            owner->b = new T[1];
            owner->c = new T[0];
        }
        public: unsigned int& operator = (const unsigned int& n) { // set function
            if(n && n != owner->_n) {
                T* a = new T[n-1];
                T* b = new T[n];
                T* c = new T[n-1];
                unsigned int n_ = n < owner->_n ? n : owner->_n;
                unsigned int _n = n_-1;
                for(unsigned int i = 0; i < _n; i++) {
                    a[i] = owner->a[i];
                    b[i] = owner->b[i];
                    c[i] = owner->c[i];
                }
                b[_n] = owner->b[_n];
                delete [] owner->a;
                delete [] owner->b;
                delete [] owner->c;
                owner->a = a;
                owner->b = b;
                owner->c = c;
                return owner->_n = n;
            } else
                return owner->_n;
        }
        public: operator unsigned int () const {                    // get function
            return owner->_n;
        }
    };
    public: _n_ n;
    public: Matrix(unsigned int n) : n(this), MatrixDiagonal<Matrix<MatrixType::TridiagonalGeneral, T>, T>(this) {
        this->n = n;
    }
    public: ~Matrix() {
        delete [] a;
        delete [] b;
        delete [] c;
    }
    public: T& operator() (const unsigned int row, const unsigned int col) {    // Matrix indexing
        if(row == col)
            return b[row];
        else if(row-1 == col)
            return a[col];
        else if(row == col-1)
            return c[row];
        other = 0;
        return other;
    }
    public: void Diagonal(const int diagonal, const T value) {
        Diagonal(diagonal, value, _n);
    }
    public: inline void Diagonal(const int diagonal, const T value, const unsigned int n) {
        unsigned int d = abs(diagonal);
        if(d <= 1)
            MatrixDiagonal<Matrix<MatrixType::TridiagonalGeneral, T>,T>::Diagonal(diagonal, value, n);
    }
    public: bool Solve(T* f, unsigned int n) {
        if(n == _n)
            return Solve(a, b, c, f, n);
        return false;
    }
    public: inline static bool Solve(T* a, T* b, T* c, T* f, unsigned int n) {     // 8n FLOPS
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
        return true;

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
template<class T> class TridiagonalElements_minus1_2_minus1 {
    public: T operator() (const unsigned int row, const unsigned int col) { // Matrix indexing
        if(row == col)
            return 2;
        else if(row-1 == col)
            return -1;
        else if(row == col-1)
            return -1;
        return 0;
    }
};
template<class T> class Matrix<MatrixType::Tridiagonal_minus1_2_minus1_6n, T> : public TridiagonalElements_minus1_2_minus1<T> {
    private: T* b;
    private: unsigned int _n;
    private: class _n_ {             // property class
        private: Matrix<MatrixType::Tridiagonal_minus1_2_minus1_6n, T>* owner;
        public: _n_(Matrix<MatrixType::Tridiagonal_minus1_2_minus1_6n, T>* owner) : owner(owner) {
            this->owner->_n = 0;
            owner->b = new T[0];
        }
        public: unsigned int& operator = (const unsigned int& n) {  // set function
            if(n && n != owner->_n) {
                delete [] owner->b;
                owner->b = new T[n];
            }
            return owner->_n = n;
        }
        public: operator unsigned int () const {                    // get function
            return owner->_n;
        }
    };
    public: _n_ n;
    public: Matrix() : n(this) {
    }
    public: Matrix(unsigned int n) : n(this) {
        this->n = n;
    }
    public: ~Matrix() {
        delete [] b;
    }
    public: static bool Solve(T* f, unsigned int n) {
        /*
         * f are source elements b_{i}, f[0] = b_{1}
         */

        T* start = f;
        T* end = &f[n-1];
        T i = 1;                        // Faster to use double than int here, and bad accuracy with int
        while(f != end)
            *f += i++/i*(*f++);         // Eq (8)
        i++;
        while(f != start) {
            *f *= 1/(i--)*i;            // Eq (11) (faster than *f /= i--/i)
            *f += *f--;                 // Eq (10)
        }
        *f /= 2;
        return true;
    }
    public: static bool Solve(T* f, unsigned int n, unsigned int cutoff) {  // 6n FLOPS => 2n FLOPS when n >> cutoff
        /*
         * f are source elements b_{i}, f[0] = b_{1}
         * The values of i >= cutoff is when (i+1)/i is approximated to 1
         */

        T* start = f;
        T* mid = &f[n < cutoff ? n-1 : cutoff-1];
        T* end = &f[n-1];
        T i = 1;                        // Faster to use double than int here, and bad accuracy with int
        while(f != mid)
            *f += i++/i*(*f++);         // Eq (8)
        while(f != end)
            *f += *f++;                 // Eq (8) with (i+1)/i = 1
        i++;
        while(f != mid)
            *f += *f--;                 // Eq (10) with (i+1)/i = 1
        while(f != start) {
            *f *= 1/(i--)*i;            // Eq (11) (faster than *f /= i--/i)
            *f += *f--;                 // Eq (10)
        }
        *f /= 2;
        return true;
    }
    public: static bool SolveInt(T* f, unsigned int n) {
        /*
         * f are source elements b_{i}, f[0] = b_{1}
         */

        T* start = f;
        T* end = &f[n-1];
        unsigned int i = 1;
        while(f != end)
            *f += (T)(i++)/i*(*f++);  // Eq (8)
        i++;
        while(f != start) {
            *f *= (T)1/(i--)*i;         // Eq (11) (faster than *f /= i--/i)
            *f += *f--;                 // Eq (10)
        }
        *f /= 2;
        return true;
    }
    public: bool SolveTrue(T* f, unsigned int n) {                       // 6n FLOPS (True)
        /*
         * f are source elements b_{i}, f[0] = b_{1}
         */

        if(_n < n)
            this->n = n;
        b[0] = 2;
        T temp;
        T* start = f;
        T* end = &f[n-1];
        while(f != end) {
            temp = (T)1/(*b++);
            *b = 2 - temp;              // Eq (9) in report
            *f += temp*(*f++);          // Eq (8)
        }
        while(f != start) {
            *f /= (*b--);               // Eq (11)
            *f += (*f--);               // Eq (10)
        }
        *f /= *b;
        return true;
    }
};
template<class T> class Matrix<MatrixType::Tridiagonal_minus1_2_minus1_4n, T> : public TridiagonalElements_minus1_2_minus1<T> {
    private: T* factor;                   // Precalculated values
    private: unsigned int _n;
    private: class _n_ {                  // property class
        private: Matrix<MatrixType::Tridiagonal_minus1_2_minus1_4n, T>* owner;
        public: _n_(Matrix<MatrixType::Tridiagonal_minus1_2_minus1_4n, T>* owner) : owner(owner) {
            this->owner->_n = 0;
            owner->factor = new T[0];
        }
        public: unsigned int & operator = (const unsigned int& n) { // set function
            if(owner->_n != n) {
                T* arr = new T[n];
                T *start;
                if(owner->_n < n) {                                 // Extend values
                    start = &arr[owner->_n]-1;
                    arr = &arr[n-1];
                    T i = n+1;
                    while(arr != start)
                        *arr-- = 1/(i--)*i;                         // Calculate (i-1)/i values
                }
                if(owner->_n != 0) {
                    start = owner->factor-1;
                    if(owner->_n < n)
                        owner->factor= &owner->factor[owner->_n-1];
                    else {
                        arr = &arr[n-1];
                        owner->factor= &owner->factor[n-1];
                    }
                    while(owner->factor != start)                   // Copy already calculated values
                        *arr-- = *(owner->factor--);
                    owner->factor++;
                }
                delete [] owner->factor;
                owner->factor = ++arr;
            }
            return owner->_n= n;
        }
        public: operator unsigned int () const {                    // get function
            return owner->_n;
        }
    };
    public: _n_ n;                              // Number of precalculated values
    public: Matrix(unsigned int n) : n(this) {
        this->n = n;
    }
    public: ~Matrix() {
        delete [] factor;
    }
    public: bool Solve(T* f, unsigned int n) {    // 4n FLOPS
        /*
         * f are source elements b_{i}, f[0] = b_{1}
         */

        if(_n < n)
            this->n = n;                // Makes enough precalculated values
        T* start = f;
        T* end = &f[n-1];
        T* factor = this->factor;
        while(f != end)
            *f += (*f++)*(*factor++);   // Eq (8)
        while(f != start) {
            *f *= *factor--;            // Eq (11)
            *f += *f--;                 // Eq (10)
        }
        *f /= 2;
        return true;
    }
};
template<class T> class Matrix<MatrixType::LU_decomposition, T> : public MatrixDiagonal<Matrix<MatrixType::LU_decomposition, T>, T> {
    public: Mat<T> matrix;
    private: unsigned int _n;
    private: class _n_ {             // property class
        private: Matrix<MatrixType::LU_decomposition, T>* owner;
        public: _n_(Matrix<MatrixType::LU_decomposition, T>* owner) : owner(owner) {
            this->owner->_n = 1;
            owner->matrix = Mat<T>(1,1, fill::zeros);
        }
        public: unsigned int& operator = (const unsigned int& n) { // set function
            if(owner->_n != n)
                owner->matrix = Mat<T>(n,n, fill::zeros);
            return owner->_n = n;
        }
        public: operator unsigned int () const {                    // get function
            return owner->_n;
        }
    };
    public: _n_ n;
    public: Matrix(unsigned int n) : n(this), MatrixDiagonal<Matrix<MatrixType::LU_decomposition, T>, T>(this) {
        this->n = n;
    }
    public: T& operator() (const unsigned int row, const unsigned int col) { // Matrix indexing
        return matrix(row, col);
    }
    public: bool Solve(T* f, unsigned int n) {
        if(n == _n) {
            Mat<T> L, U;
            lu(L, U, matrix);                   // Armadillo LU decomposition
            int i = 0, j, n = matrix.n_rows;
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
            return true;
        }
        return false;
    }
};
#endif // MATRIX_H
