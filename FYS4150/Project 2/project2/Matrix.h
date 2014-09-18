/*
 *  FYS4150 - Computational Physics - Project 1
 *
 *  Written by: Eimund Smestad
 *
 *  16.09.2014
 *
 *  c11 compiler
 */

#ifndef MATRIX_H
#define MATRIX_H

#include "armadillo"
#include "lib.cpp"

using namespace arma;

enum class MatrixType {
    Symmetric,
    Tridiagonal,
    Tridiagonal_m1_2_m1,
    Tridiagonal_m1_2_m1_6n,
    Tridiagonal_m1_2_m1_4n,
    LU_decomposition
};

template<MatrixType Type, class T> class Matrix;
template<class M, class T> class MatrixDiagonal {
    private: M* owner;
    private: template<class P> class Type {
        public: static inline void Diagonal(M* owner, const int diagonal, const P value, const unsigned int n) {
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
    private: template<class P> class Type<P*> {
        public: static inline void Diagonal(M* owner, const int diagonal, P* value, const unsigned int n) {
            unsigned int d = abs(diagonal);
            owner->n = n;
            unsigned int nmax = n - d;
            if(diagonal > 0) {
                for(unsigned int i = 0; i < nmax; i++)
                    owner->operator()(i,i+d) = value[i];
            } else {
                for(unsigned int i = 0; i < nmax; i++)
                    owner->operator()(i+d,i) = value[i];
            }
        }
    };
    protected: MatrixDiagonal(M* owner) : owner(owner) {
    }
    public: template<class P> inline void Diagonal(const int diagonal, const P value) {
        Type<P>::Diagonal(owner, diagonal, value, owner->n);
    }
    public: template<class P> inline void Diagonal(const int diagonal, const P value, const unsigned int n) {
        Type<P>::Diagonal(owner, diagonal, value, n);
    }
};
template <class M, class T> class MatrixElements {
    private: T** matrix;
    private: unsigned int _n;
    private: class _n_ {             // property class
        private: M* owner;
        public: _n_(M* owner) : owner(owner) {
            this->owner->_n = 1;
            owner->matrix = new T*[1];
            owner->matrix[0] = new T[1];
        }
        public: unsigned int& operator = (const unsigned int& n) { // set function
            if(n && n != owner->_n) {
                T** matrix = new T*[n];
                for(unsigned int i = 0; i < n; i++)
                    matrix[i] = new T[n];
                unsigned int n_ = n < owner->_n ? n : owner->_n;
                for(unsigned int i = 0; i < n_; i++) {
                    for(unsigned int j = 0; j < n; j++)
                        matrix[i][j] = owner->matrix[i][j];
                }
                for(unsigned int i = 0; owner->_n; i++)
                    delete [] matrix[i];
                delete [] matrix;
                owner->matrix = matrix;
                return owner->_n = n;
            } else
                return owner->_n;
        }
        public: operator unsigned int () const {                    // get function
            return owner->_n;
        }
    };
    public: _n_ n;
    protected: MatrixElements(unsigned int n) : n(this) {
        this->n = n;
    }
    public: ~MatrixElements() {
        for(unsigned int i = 0; i < _n; i++)
            delete [] matrix[i];
        delete [] matrix;
    }
    public: T& operator() (const unsigned int row, const unsigned int col) { // Matrix indexing
        return matrix[row][col];
    }
};
template <class T> class MatrixElements<Matrix<MatrixType::Symmetric, T>, T> {
    private: T** matrix;
    private: unsigned int _n;
    private: class _n_ {             // property class
        private: MatrixElements<Matrix<MatrixType::Symmetric, T>, T>* owner;
        public: _n_(MatrixElements<Matrix<MatrixType::Symmetric, T>, T>* owner) : owner(owner) {
            this->owner->_n = 1;
            owner->matrix = new T*[1];
            owner->matrix[0] = new T[1];
        }
        public: unsigned int& operator = (const unsigned int& n) { // set function
            if(n && n != owner->_n) {
                T** matrix = new T*[n];
                for(unsigned int i = 0; i < n; i++)
                    matrix[i] = new T[n-i];
                unsigned int n_ = n < owner->_n ? n : owner->_n;
                for(unsigned int i = 0; i < n_; i++) {
                    for(unsigned int j = 0; j < n_-i; j++)
                        matrix[i][j] = owner->matrix[i][j];
                }
                for(unsigned int i = 0; i < owner->_n; i++)
                    delete [] owner->matrix[i];
                delete [] owner->matrix;
                owner->matrix = matrix;
                return owner->_n = n;
            } else
                return owner->_n;
        }
        public: operator unsigned int () const {                    // get function
            return owner->_n;
        }
    };
    public: _n_ n;
    protected: MatrixElements(unsigned int n) : n(this) {
        this->n = n;
    }
    public: ~MatrixElements() {
        for(unsigned int i = 0; i < _n; i++)
            delete [] matrix[i];
        delete [] matrix;
    }
    public: T& operator() (const unsigned int row, const unsigned int col) { // Matrix indexing
        return row < col ? matrix[row][col-row] : matrix[col][row-col];
    }
};
template <class T> class MatrixElements<Matrix<MatrixType::Tridiagonal, T>, T> {
    private: T other;
    protected: T *a, *b, *c;
    protected: unsigned int _n;
    private: class _n_ {             // property class
        private: MatrixElements<Matrix<MatrixType::Tridiagonal, T>, T>* owner;
        public: _n_(MatrixElements<Matrix<MatrixType::Tridiagonal, T>, T>* owner) : owner(owner) {
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
    public: MatrixElements(unsigned int n) : n(this) {
        this->n = n;
    }
    public: ~MatrixElements() {
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
};
template <class T> class MatrixElements<Matrix<MatrixType::Tridiagonal_m1_2_m1, T>, T> {
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
template<class T> class Matrix<MatrixType::Symmetric, T> :
        public MatrixDiagonal<Matrix<MatrixType::Symmetric, T>, T>,
        public MatrixElements<Matrix<MatrixType::Symmetric, T>, T> {
    public: Matrix(unsigned int n) :
            MatrixDiagonal<Matrix<MatrixType::Symmetric, T>, T>(this),
            MatrixElements<Matrix<MatrixType::Symmetric, T>, T>(n) {
    }
};
template<class T> class Matrix<MatrixType::Tridiagonal, T> :
        public MatrixDiagonal<Matrix<MatrixType::Tridiagonal, T>, T>,
        public MatrixElements<Matrix<MatrixType::Tridiagonal, T>, T> {
    public: Matrix(unsigned int n) :
        MatrixDiagonal<Matrix<MatrixType::Tridiagonal, T>, T>(this),
        MatrixElements<Matrix<MatrixType::Tridiagonal, T>, T>(n) {
    }
    public: ~Matrix() {
    }
    public: void Diagonal(const int diagonal, const T value) {
        Diagonal(diagonal, value, MatrixElements<Matrix<MatrixType::Tridiagonal, T>, T>::_n);
    }
    public: inline void Diagonal(const int diagonal, const T value, const unsigned int n) {
        unsigned int d = abs(diagonal);
        if(d <= 1)
            MatrixDiagonal<Matrix<MatrixType::Tridiagonal, T>,T>::Diagonal(diagonal, value, n);
    }
    public: bool Solve(T* f, unsigned int n) {
        if(n == MatrixElements<Matrix<MatrixType::Tridiagonal, T>, T>::_n)
            return Solve(MatrixElements<Matrix<MatrixType::Tridiagonal, T>, T>::a, MatrixElements<Matrix<MatrixType::Tridiagonal, T>, T>::b, MatrixElements<Matrix<MatrixType::Tridiagonal, T>, T>::c, f, n);
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
            *b -= temp*(*c++);          // Eq (6) in report
            *f -= temp*(*f++);          // Eq (5)
        }
        while(f != start) {
            *f /= (*b--);               // Eq (9)
            *f -= (*--c)*(*f--);        // Eq (8)
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
template<class T> class Matrix<MatrixType::Tridiagonal_m1_2_m1_6n, T> :
        public MatrixElements<Matrix<MatrixType::Tridiagonal_m1_2_m1, T>, T>  {
    public: static bool Solve(T* f, unsigned int n) {
        /*
         * f are source elements b_{i}, f[0] = b_{1}
         */

        T* start = f;
        T* end = &f[n-1];
        T i = 1;                        // Faster to use double than int here, and bad accuracy with int
        while(f != end)
            *f += i++/i*(*f++);         // Eq (10)
        i++;
        while(f != start) {
            *f *= 1/(i--)*i;            // Eq (13) (faster than *f /= i--/i)
            *f += *f--;                 // Eq (12)
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
            *f += i++/i*(*f++);         // Eq (10)
        while(f != end)
            *f += *f++;                 // Eq (10) with (i+1)/i = 1
        i++;
        while(f != mid)
            *f += *f--;                 // Eq (12) with (i+1)/i = 1
        while(f != start) {
            *f *= 1/(i--)*i;            // Eq (13) (faster than *f /= i--/i)
            *f += *f--;                 // Eq (12)
        }
        *f /= 2;
        return true;
    }
};
template<class T> class Matrix<MatrixType::Tridiagonal_m1_2_m1_4n, T> :
        public MatrixElements<Matrix<MatrixType::Tridiagonal_m1_2_m1, T>, T>  {
    private: T* factor;                   // Precalculated values
    private: unsigned int _n;
    private: class _n_ {                  // property class
        private: Matrix<MatrixType::Tridiagonal_m1_2_m1_4n, T>* owner;
        public: _n_(Matrix<MatrixType::Tridiagonal_m1_2_m1_4n, T>* owner) : owner(owner) {
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
            *f += (*f++)*(*factor++);   // Eq (10)
        while(f != start) {
            *f *= *factor--;            // Eq (13)
            *f += *f--;                 // Eq (12)
        }
        *f /= 2;
        return true;
    }
};
template<class T> class Matrix<MatrixType::LU_decomposition, T> :
        public MatrixDiagonal<Matrix<MatrixType::LU_decomposition, T>, T> {
    public: Mat<T> matrix, L, U;
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
    public: void LU() {
        lu(L, U, matrix);                   // Armadillo LU decomposition
    }
    public: bool Solve(T* f, unsigned int n) {
        if(n == _n) {
            T* f_;
            int i = 0, j;
            while(++i < n) {                // Forward solve Ly = f
                f_ = &f[i];
                for(j = 0; j < i; j++)
                    *f_ -= L(i,j)*f[j];
            }
            n--;
            while(i--) {                    // Backward solve Ux = y
                f_ = &f[i];
                for(j = n; j > i; j--)
                    *f_ -= U(i,j)*f[j];
                *f_ /= U(i,i);
            }
            return true;
        }
        return false;
    }
};
#endif // MATRIX_H
