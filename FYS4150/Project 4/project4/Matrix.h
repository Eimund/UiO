/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  05.11.2014
 *
 *  c11 compiler
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <iostream>
#include "armadillo"
#include "lib.cpp"
#include "Array.h"
#include "Delegate.h"
#include "Property.h"

using namespace arma;
using namespace std;

enum class MatrixType {
    Square,
    SquareT,
    Symmetric,
    Tridiagonal,
    TridiagonalSymmetric,
    Tridiagonal_m1_X_m1,
    Tridiagonal_m1_C_m1,
    Tridiagonal_m1_2_m1,
    Tridiagonal_m1_2_m1_6n,
    Tridiagonal_m1_2_m1_4n,
    LU_decomposition,
    tqli
};
typedef struct {
    unsigned int i, j;
} MatrixIndex;

template<MatrixType Type, class T> class Matrix;
template<MatrixType Type, class T> void MatrixCout(Matrix<Type, T>& matrix) {
    for(unsigned int i = 0; i < matrix.n; i++) {
        for(unsigned int j = 0; j < matrix.n; j++)
            cout << matrix(i,j) << "\t";
        cout << '\n';
    }
    cout << '\n';
}
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
        public: static inline void Diagonal(M* owner, const int diagonal, const P* value, const unsigned int n) {
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
    protected: T** matrix;
    protected: unsigned int _n;
    private: class _n_ {             // property class
        private: MatrixElements<M, T>* owner;
        public: _n_(MatrixElements<M, T>* owner) : owner(owner) {
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
                    for(unsigned int j = 0; j < n_; j++)
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
    public: inline T& operator() (const unsigned int& row, const unsigned int& col) { // Matrix indexing
        return matrix[row][col];
    }
    public: inline operator T**() const {
        return matrix;
    }
    public: void Clear() {
        for(unsigned int i = 0; i < _n; i++) {
            for(unsigned int j = 0; j < _n; j++)
                matrix[i][j] = 0;
        }
    }
};
template <class T> class MatrixElements<Matrix<MatrixType::SquareT, T>, T> {
    protected: T** matrix;
    protected: unsigned int _n;
    private: class _n_ {             // property class
        private: MatrixElements<Matrix<MatrixType::SquareT, T>, T>* owner;
        public: _n_(MatrixElements<Matrix<MatrixType::SquareT, T>, T>* owner) : owner(owner) {
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
                    for(unsigned int j = 0; j < n_; j++)
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
    public: inline T& operator() (const unsigned int& row, const unsigned int& col) { // Matrix indexing
        return matrix[col][row];
    }
    public: inline operator T**() const {
        return matrix;
    }
    public: void Clear() {
        for(unsigned int i = 0; i < _n; i++) {
            for(unsigned int j = 0; j < _n; j++)
                matrix[i][j] = 0;
        }
    }
    public: void Normailze() {
        T sum;
        for(unsigned int i = 0, j; i < n; i++) {
            sum = 0;
            for(j = 0; j < n; j++)
                sum += matrix[i][j];
            sum = (T)1/sqrt(sum);
            for(j = 0; j < n; j++)
                matrix[i][j] *= sum;
        }
    }
    public: void Transpose() {
        T temp;
        for(unsigned int i = 0; i < n; i++) {
            for(unsigned int j = i; j < n; j++) {
                temp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = temp;
            }
        }
    }
};
template <class T> class MatrixElements<Matrix<MatrixType::Symmetric, T>, T> {
    protected: T** matrix;
    protected: unsigned int _n;
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
    public: inline T& operator() (const unsigned int& row, const unsigned int& col) { // Matrix indexing
        return row < col ? matrix[col-row][row] : matrix[row-col][col];
    }
    public: inline operator T**() const {
        return matrix;
    }
    public: void Clear() {
        for(unsigned int i = 0, n = _n; i < _n; i++, n--) {
            for(unsigned int j = 0; j < n; j++)
                matrix[i][j] = 0;
        }
    }
};
template <class T> class MatrixElements<Matrix<MatrixType::tqli, T>, T> {
    protected: T** matrix;
    protected: unsigned int _n;
    private: class _n_ {             // property class
        private: MatrixElements<Matrix<MatrixType::tqli, T>, T>* owner;
        public: _n_(MatrixElements<Matrix<MatrixType::tqli, T>, T>* owner) : owner(owner) {
            this->owner->_n = 1;
            owner->matrix = new T*[1];
            owner->matrix[0] = new T[1];
        }
        public: unsigned int& operator = (const unsigned int& n) { // set function
            if(n && n != owner->_n) {
                T** matrix = new T*[n];
                for(unsigned int i = 0; i < n; i++) {
                    if(i == 1) {
                        matrix[1] = new T[n];           // Bug i tqli i lib.cpp, trenger ett ekstra element
                        matrix[1] = &matrix[1][1];
                    }
                    else
                        matrix[i] = new T[n-i];
                }
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
        for(unsigned int i = 0; i < _n; i++) {
            if(i == 1)
                delete [] &matrix[1][-1];           // Bug i tqli i lib.cpp,
            else
                delete [] matrix[i];
        }
        delete [] matrix;
    }
    public: inline T& operator() (const unsigned int& row, const unsigned int& col) { // Matrix indexing
        return row < col ? matrix[col-row][row] : matrix[row-col][col];
    }
    public: void Clear() {
        for(unsigned int i = 0, n = _n; i < _n; i++, n--) {
            for(unsigned int j = 0; j < n; j++)
                matrix[i][j] = 0;
        }
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
template <class T> class MatrixElements<Matrix<MatrixType::TridiagonalSymmetric, T>, T> {
    private: T other;
    protected: T *a, *b;
    protected: unsigned int _n;
    private: class _n_ {             // property class
        private: MatrixElements<Matrix<MatrixType::TridiagonalSymmetric, T>, T>* owner;
        public: _n_(MatrixElements<Matrix<MatrixType::TridiagonalSymmetric, T>, T>* owner) : owner(owner) {
            this->owner->_n = 1;
            owner->a = new T[0];
            owner->b = new T[1];
        }
        public: unsigned int& operator = (const unsigned int& n) { // set function
            if(n && n != owner->_n) {
                T* a = new T[n-1];
                T* b = new T[n];
                unsigned int n_ = n < owner->_n ? n : owner->_n;
                unsigned int _n = n_-1;
                for(unsigned int i = 0; i < _n; i++) {
                    a[i] = owner->a[i];
                    b[i] = owner->b[i];
                }
                b[_n] = owner->b[_n];
                delete [] owner->a;
                delete [] owner->b;
                owner->a = a;
                owner->b = b;
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
    }
    public: T& operator() (const unsigned int row, const unsigned int col) {    // Matrix indexing
        if(row == col)
            return b[row];
        else if(row-1 == col)
            return a[col];
        else if(row == col-1)
            return a[row];
        other = 0;
        return other;
    }
};
template <class T> class MatrixElements<Matrix<MatrixType::Tridiagonal_m1_X_m1, T>, T> {
    private: T a, other;
    protected: T *b;
    protected: unsigned int _n;
    private: class _n_ {             // property class
        private: MatrixElements<Matrix<MatrixType::Tridiagonal_m1_X_m1, T>, T>* owner;
        public: _n_(MatrixElements<Matrix<MatrixType::Tridiagonal_m1_X_m1, T>, T>* owner) : owner(owner) {
            this->owner->_n = 1;
            owner->b = new T[1];
        }
        public: unsigned int& operator = (const unsigned int& n) { // set function
            if(n && n != owner->_n) {
                T* b = new T[n];
                unsigned int n_ = n < owner->_n ? n : owner->_n;
                for(unsigned int i = 0; i < n_; i++)
                    b[i] = owner->b[i];
                delete [] owner->b;
                owner->b = b;
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
        delete [] b;
    }
    public: T& operator() (const unsigned int row, const unsigned int col) {    // Matrix indexing
        a = -1;
        if(row == col)
            return b[row];
        else if(row-1 == col)
            return a;
        else if(row == col-1)
            return a;
        other = 0;
        return other;
    }
};
template <class T> class MatrixElements<Matrix<MatrixType::Tridiagonal_m1_C_m1, T>, T> {
    private: T a, other;
    protected: T b;
    public: unsigned int n;
    public: MatrixElements() = default;
    public: MatrixElements(unsigned int n) : n(n) {
    }
    public: T& operator() (const unsigned int row, const unsigned int col) { // Matrix indexing
        a = -1;
        if(row == col)
            return b;
        else if(row-1 == col)
            return a;
        else if(row == col-1)
            return a;
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
template<class T> class Matrix<MatrixType::Square, T> :
        public MatrixDiagonal<Matrix<MatrixType::Square, T>, T>,
        public MatrixElements<Matrix<MatrixType::Square, T>, T> {
    public: Matrix(unsigned int n) :
        MatrixDiagonal<Matrix<MatrixType::Square, T>, T>(this),
        MatrixElements<Matrix<MatrixType::Square, T>, T>(n) {

    }
	void JacobiMethod(T error, unsigned int &num /**/) {
        MatrixIndex I;
        T s, c, t, tau, a_kk, a_ll, a_ik, a_il;
        num = 0;                            // To be removed
        while(error < MaxAbs(this->matrix, this->_n, I)) {
            tau = (this->matrix[I.j][I.j] - this->matrix[I.i][I.i])/(2*this->matrix[I.i][I.j]);
            if (tau > 0)
                t = 1.0/(tau + sqrt(1.0 + tau*tau));
			else
                t = -1.0/( -tau + sqrt(1.0 + tau*tau));
            c = 1/sqrt(1+t*t);
            s = c*t;

            a_kk = this->matrix[I.i][I.i];
            a_ll = this->matrix[I.j][I.j];
            // changing the matrix elements with indices k and l
            this->matrix[I.i][I.i] = c*c*a_kk - 2.0*c*s*this->matrix[I.i][I.j] + s*s*a_ll;
            this->matrix[I.j][I.j] = s*s*a_kk + 2.0*c*s*this->matrix[I.i][I.j] + c*c*a_ll;
            this->matrix[I.i][I.j] = 0.0; // hard-coding of the zeros
            this->matrix[I.j][I.i] = 0.0;
            // and then we change the remaining elements
            for (int i = 0; i < this->_n; i++ ) {
                if ( i != I.i && i != I.j ) {
                    a_ik = this->matrix[i][I.i];
                    a_il = this->matrix[i][I.j];
                    this->matrix[i][I.i] = c*a_ik + s*a_il;
                    this->matrix[I.i][i] = this->matrix[i][I.i];
                    this->matrix[i][I.j] = c*a_il - s*a_ik;
                    this->matrix[I.j][i] = this->matrix[i][I.j];
				}
            }
            num++;                      // To be removed
		}
	}
	public: T MaxAbs(T** matrix, unsigned int n, MatrixIndex& I) {
		T max = 0, abs1;
		for(unsigned int i = 0, j; i < n; i++) {
			for(j = i+1; j < n; j++) {
				abs1 = abs(matrix[i][j]);
				if(abs1 > max) {
					max = abs1;
					I.i = i;
					I.j = j;
				}
			}
		}
		return max;
	}
};
template<class T> class Matrix<MatrixType::SquareT, T> :
        public MatrixDiagonal<Matrix<MatrixType::SquareT, T>, T>,
        public MatrixElements<Matrix<MatrixType::SquareT, T>, T> {
    public: Matrix(unsigned int n) :
        MatrixDiagonal<Matrix<MatrixType::SquareT, T>, T>(this),
        MatrixElements<Matrix<MatrixType::SquareT, T>, T>(n) {
    }
};
template<class T> class Matrix<MatrixType::Symmetric, T> :
		public MatrixDiagonal<Matrix<MatrixType::Symmetric, T>, T>,
        public MatrixElements<Matrix<MatrixType::Symmetric, T>, T> {
    public: Matrix(unsigned int n) :
            MatrixDiagonal<Matrix<MatrixType::Symmetric, T>, T>(this),
            MatrixElements<Matrix<MatrixType::Symmetric, T>, T>(n) {
	}
	public: T* JacobiMethod(T error, unsigned int& num /**/) {
		MatrixIndex I;
		T** offdiagonal = &this->matrix[1];
		unsigned int offn = this->_n-1;
		num = 0;                                    // To be removed
		while(error < MaxAbs(offdiagonal, offn, I)) // Find the maximum element, and stop when error is smaller than desired
			Rotate(I.j, ++I.i, num);
		return this->matrix[0];
	}
	public: T* JacobiMethodFD(T error, unsigned int& num /**/) {        // Forward diagonal iteration
		bool run = true;
		num = 0;                                                        // To be removed
		while(run) {
			run = false;
			for(unsigned int d = 1, r, n = this->_n-1; n; d++, n--) {   // Run through the diagonals
				for(r = 0; r < n; r++) {                                // Run through the rows
					if(error < abs(this->matrix[d][r])) {
						run = true;
						Rotate(r, d, num);
					}
				}
			}
		}
		return this->matrix[0];
	}
	public: T* JacobiMethodRD(T error, unsigned int& num /**/) {        // Reverse diagonal iteration
		bool run = true;
		num = 0;                                                        // To be removed
		while(run) {
			run = false;
			for(unsigned int d = this->_n-1, r, p = 1; d; d--, p++) {   // Run through the diagonals
				for(r = 0; r < p; r++) {                                // Run through the rows
					if(error < abs(this->matrix[d][r])) {
						run = true;
						Rotate(r, d, num);
					}
				}
			}
		}
		return this->matrix[0];
	}
	public: T* JacobiMethodFC(T error, unsigned int& num /**/) {        // Forward column iteration
		bool run;
		num = 0;                                                        // To be removed
		do {
			run = false;
			for(unsigned int r = 0, d, n = this->_n; r < this->_n; r++, n--) {
				for(d = 1; d < n; d++) {
					if(error < abs(this->matrix[d][r])) {
						run = true;
						Rotate(r, d, num);
					}
				}
			}
		} while(run);
		return this->matrix[0];
	}
	public: T* JacobiMethodRC(T error, unsigned int& num /**/) {        // Reverse column iteration
		bool run;
		num = 0;                                                        // To be removed
		do {
			run = false;
			for(unsigned int r = 0, d, n = this->_n-1; r < this->_n; r++, n--) {
				for(d = n; d; d--) {
					if(error < abs(this->matrix[d][r])) {
						run = true;
						Rotate(r, d, num);
					}
				}
			}
		} while(run);
		return this->matrix[0];
	}
	public: T* JacobiMethodFR(T error, unsigned int& num /**/) {        // Forward row iteration
		bool run;
		num = 0;                                                        // To be removed
		do {
			run = false;
			for(unsigned int r, d = 1, c; d < this->_n; d++) {
				for(r = 0, c = d; c; r++, c--) {
					if(error < abs(this->matrix[c][r])) {
						run = true;
						Rotate(r, c, num);
					}
				}
			}
		} while(run);
		return this->matrix[0];
	}
	public: T* JacobiMethodRR(T error, unsigned int& num /**/) {        // Reverse row iteration
		bool run;
		num = 0;                                                        // To be removed
		do {
			run = false;
			for(unsigned int r, d = 1, c; d < this->_n; d++) {
				for(r = d-1, c = 1; c <= d; r--, c++) {
					if(error < abs(this->matrix[c][r])) {
						run = true;
						Rotate(r, c, num);
					}
				}
			}
		} while(run);
		return this->matrix[0];
	}
	private: T MaxAbs(T** matrix, unsigned int n, MatrixIndex& I) {
		T max = 0, abs1;
		for(unsigned int i = 0, j; n; i++, n--) {
			for(j = 0; j < n; j++) {
				abs1 = abs(matrix[i][j]);
				if(abs1 > max) {
					max = abs1;
					I.i = i;            // I.i is the diagonal number = col-row
					I.j = j;            // I.j is the row number
				}
			}
		}
		return max;
	}
	private: void Rotate(unsigned int i, unsigned int j, unsigned int& num) {
        T b = 2*this->matrix[j][i];
        this->matrix[j][i] = 0;     // bij=bij=0

		j += i;                     // Convert j from diagonalnumber to column number

        T a = this->matrix[0][j] - this->matrix[0][i];  // Calculate cos and sin
        T t = a/b;
        if(t > 0)
			t = (T)1/(t+sqrt(1+t*t));   // Calculate -cot !
		else
			t = (T)1/(t-sqrt(1+t*t));   // Calculate -cot !
		T s2 = (T)1/(1+t*t);
        T s = sqrt(s2);
        T c = t*s;

		a = a*s2 + b*s*c;               // s*c = - sin*cos because of t=-cot
		this->matrix[0][i] += a;       // bii eq(18)
		this->matrix[0][j] -= a;       // bjj eq(19)

        unsigned int k, l, n;
        for(k = i, l = j, n = 0; n < i; k--, l--, n++) {    // Calculate bik and bjk, row < i and row < j
            a = this->matrix[k][n];
            b = this->matrix[l][n];
            this->matrix[k][n] = a*c - b*s;
            this->matrix[l][n] = b*c + a*s;
         }
		 for(k++, l--, n++; n < j; k++, l--, n++) {         // row > i and row < j
            a = this->matrix[k][i];
            b = this->matrix[l][n];
            this->matrix[k][i] = a*c - b*s;
            this->matrix[l][n] = b*c + a*s;
        }
        for(k++, l++, n = this->_n-j; l < n; k++, l++) {    // row > i and row > j
            a = this->matrix[k][i];
			b = this->matrix[l][j];
            this->matrix[k][i] = a*c - b*s;
            this->matrix[l][j] = b*c + a*s;
        }
        num++;                                              // To be removed
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
    }
};
template<class T> class Matrix<MatrixType::TridiagonalSymmetric, T> :
		public MatrixDiagonal<Matrix<MatrixType::TridiagonalSymmetric, T>, T>,
        public MatrixElements<Matrix<MatrixType::TridiagonalSymmetric, T>, T> {
    public: Matrix(unsigned int n) :
        MatrixDiagonal<Matrix<MatrixType::TridiagonalSymmetric, T>, T>(this),
        MatrixElements<Matrix<MatrixType::TridiagonalSymmetric, T>, T>(n) {
    }
    public: ~Matrix() {
    }
    public: template<class P> void Diagonal(const int diagonal, const P value) {
        Diagonal(diagonal, value, MatrixElements<Matrix<MatrixType::TridiagonalSymmetric, T>, T>::_n);
    }
	public: template<class P> inline void Diagonal(const int diagonal, const P value, const unsigned int n) {
        unsigned int d = abs(diagonal);
		if(d <= 1)
            MatrixDiagonal<Matrix<MatrixType::TridiagonalSymmetric, T>, T>::Diagonal(diagonal, value, n);
	}
	public: T* QRalgorithm(T error, unsigned int& num) {
        unsigned int i, j, n = this->_n-2;
        T c, s, _c, _s, temp, temp2, max, a, b;

        num = 0;
		do {
			// First rotation in QR transformation
            max = 0;
            a = this->a[0];
            b = this->b[0];

            temp = sqrt((T)1/(a*a+b*b));            // Calculate sin and cos
            c = b*temp;
            s = a*temp;

            temp = a*s;                             // Eigenvalues
            this->b[0] = b*c*c + this->b[1]*s*s + 2*temp*c;
            this->b[1] = this->b[1]*c - temp;

            num++;

			// The intermediate rotations in QR transformation
            i = 1;
            while(i < n) {
                if(abs(this->a[i]) > error) {       // Offdiagonal element still to large?

                    _c = c;
                    _s = s;
                    a = this->a[i];
                    b = this->b[i];

                    temp = (T)1/sqrt(a*a + b*b);    // Calculate sin and cos
                    c = b*temp;
                    s = a*temp;

                    temp = a*s;                     // Eigenvalues
                    b = b*c + temp;
                    this->a[i-1] = b*_s;
                    temp2 = abs(this->a[i-1]);
                    if(max < temp2)
                        max = temp2;
                    temp *= _c;
                    this->b[i] = b*_c*c + temp*c + this->b[++i]*s*s;
                    this->b[i] = this->b[i]*c - temp;

                    num++;

                } else {    // Offdiagonal element less then error
                    if(s) { // If previous sin not zero, then update elements
                        _c = c;
                        _s = s;

                        this->a[i-1] = this->b[i]*_s;   // Eigenvalues
                        this->b[i] = this->b[i]*_c;
                        temp2 = abs(this->a[i-1]);
                        if(max < temp2)
                            max = temp2;

                        c = 1;
                        s = 0;
                    }
					i++;
                }
			}

			// Last rotation in QR transformation
            a = this->a[i];
            b = this->b[i];

            temp = (T)1/sqrt(a*a + b*b);    // Calculate sin and cos
            _c = b*temp;
            _s = a*temp;

            temp = a*_s;                    // Eigenvalues
            b = b*_c + temp;
			this->a[i-1] = b*s;
            temp2 = abs(this->a[i-1]);
            if(max < temp2)
                max = temp2;
            temp *= c;
            this->b[i] = b*_c*c + temp*_c + this->b[++i]*_s*_s;
            b = this->b[i]*_c - temp;
            this->a[i-1] = b*_s;
            temp2 = abs(this->a[i-1]);
            if(max < temp2)
                max = temp2;
            this->b[i] = b*_c;

            num++;

        } while(max > error);

        return this->b;
    }
};
template<class T> class Matrix<MatrixType::Tridiagonal_m1_X_m1, T> :
        public MatrixDiagonal<Matrix<MatrixType::Tridiagonal_m1_X_m1, T>, T>,
        public MatrixElements<Matrix<MatrixType::Tridiagonal_m1_X_m1, T>, T>  {
    public: Matrix(unsigned int n) :
        MatrixDiagonal<Matrix<MatrixType::Tridiagonal_m1_X_m1, T>, T>(this),
        MatrixElements<Matrix<MatrixType::Tridiagonal_m1_X_m1, T>, T>(n) {
    }
    public: template<class P> void Diagonal(const P value) {
        Diagonal(0, value, MatrixElements<Matrix<MatrixType::Tridiagonal_m1_X_m1, T>, T>::_n);
    }
    public: template<class P> void Diagonal(const int diagonal, const P value) {
        Diagonal(diagonal, value, MatrixElements<Matrix<MatrixType::Tridiagonal_m1_X_m1, T>, T>::_n);
    }
    public: template<class P> inline void Diagonal(const int diagonal, const P value, const unsigned int n) {
        if(diagonal == 0)
            MatrixDiagonal<Matrix<MatrixType::Tridiagonal_m1_X_m1, T>, T>::Diagonal(diagonal, value, n);
    }
    public: T* Eigenvector(T value) {
        T temp = 1, temp2 = this->b[0] - value;
        this->b[0] = 1;
        for(unsigned int i = 1; i < this->_n; i++) {
            temp  *= temp2;
            temp2 = this->b[i] - value - (T)1/temp2;
            this->b[i] = temp;
        }
        return this->b;
    }
    public: bool Solve(T* f, unsigned int n) {
        if(this->_n < n)
            this->n = n;
        T temp;
        unsigned int i, i1;
        for(i = 0, i1 = 1; i1 < n; i++, i1++) {
            temp = (T)1/this->b[i];
            this->b[i1] -= temp;
            f[i1] += temp*f[i];
        }
        for(i--, i1--; i1; i--, i1--) {
            f[i1] /= this->b[i1];
            f[i] += f[i1];
        }
        f[0] /= this->b[0];
        return true;
    }
    public: bool Solve(T** f, unsigned int n1, unsigned int n2) {
        if(this->_n < n2)
            this->n = n2;
        T temp;
        unsigned int i, i1, j;
        for(i = 0, i1 = 1; i1 < n2; i++, i1++) {
            temp = (T)1/this->b[i];
            this->b[i1] -= temp;
            for(j = 0; j < n1; j++)
                f[j][i1] += temp*f[j][i];
        }
        for(i--, i1--; i1; i--, i1--) {
            for(j = 0; j < n1; j++) {
                f[j][i1] /= this->b[i1];
                f[j][i] += f[j][i1];
            }
        }
        for(j = 0; j < n1; j++)
            f[j][0] /= this->b[0];
        return true;
    }
};
template<class T> class Matrix<MatrixType::Tridiagonal_m1_C_m1, T> :
        public MatrixElements<Matrix<MatrixType::Tridiagonal_m1_C_m1, T>, T>  {
    private: typedef Matrix<MatrixType::Tridiagonal_m1_C_m1, T> THIS;
    private: T* factor;                   // Precalculated values
    public: Property<PropertyType::ReadOnly, THIS, ArrayLength<THIS, T>, int> n;
    public: Matrix() = default;
    public: Matrix(unsigned int n) :
        MatrixElements<THIS, T>(n),
        n(PROPERTY(this->n, factor, n, Delegate<THIS, void, T*, unsigned int, unsigned int>(this, &THIS::SolveInitialize))) {
    }
    public: ~Matrix() {
        delete [] factor;
    }
    public: template<class P> void Diagonal(const P value) {
        Diagonal(0, value, 0);
    }
    public: template<class P> void Diagonal(const int diagonal, const P value) {
        Diagonal(diagonal, value, 0);
    }
    public: template<class P> inline void Diagonal(const int diagonal, const P value, unsigned int) {
        if(diagonal == 0) {
            this->b = value;
            SolveInitialize();
        }
    }
    public: bool Solve(T* f, int n) {
        if(this->n < n)
            this->n = n;
        T* start = f;
        T* end = &f[n-1];
        T* factor = this->factor;
        while(f != end)
            *f += (*f++)*(*factor++);
        while(f != start) {
            *f *= *factor--;
            *f += *f--;
        }
        *f /= 2;
        return true;
    }
    public: void SolveInitialize() {
        SolveInitialize(factor, 0, n);
    }
    private: inline void SolveInitialize(T* array, unsigned int i, unsigned int n) {
        MatrixElements<Matrix<MatrixType::Tridiagonal_m1_C_m1, T>, T>:: n = n;
        if(i == 0 && n > 0) {
            array[0] = this->b;
            i++;
        }
        while(array[i-1] && i < n) {
            array[i] = this->b - (T)1/array[i-1];
            i++;
        }
    }
    public: inline void SolveWithInitialize(T* f, unsigned int n) {
        SolveInitialize(factor, 0, n);
        Solve(f, n);
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
template<class T> class Matrix<MatrixType::tqli, T> :
        public MatrixDiagonal<Matrix<MatrixType::tqli, T>, T>,
        public MatrixElements<Matrix<MatrixType::tqli, T>, T> {
    public: Matrix(unsigned int n) :
        MatrixDiagonal<Matrix<MatrixType::tqli, T>, T>(this),
        MatrixElements<Matrix<MatrixType::tqli, T>, T>(n) {
    }
    public: T* tqli(T** z) {
        ::tqli(this->matrix[0], &this->matrix[1][-1], this->_n, z);
        return this->matrix[0];
    }
};
#endif // MATRIX_H
