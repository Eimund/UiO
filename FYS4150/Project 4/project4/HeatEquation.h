/*
 *  FYS4150 - Computational Physics - Project 4
 *
 *  Written by: Eimund Smestad
 *
 *  07.11.2014
 *
 *  c11 compiler
 */
#ifndef HEATEQUATION_H
#define HEATEQUATION_H

#include <fstream>
#include "Array.h"
#include "Boundary.h"
#include "Matrix.h"

template<typename T, unsigned int D> class HeatEquation : public Boundary<T,D> {
    private: typedef HeatEquation<T, D> THIS;
    private: T _1_theta;
    private: T alpha[D];
    public: T* u[D];
    public: Matrix<MatrixType::Tridiagonal_m1_C_m1, T>* matrix[D];
    public: Property<PropertyType::ReadOnly, THIS, PropertyGet<THIS, ArrayLength<Boundary<T,D>,T>>, ArrayLength<Boundary<T,D>,T>> X[D+1];
    public: Property<PropertyType::ReadOnly, THIS, PropertySetGet<THIS, int>, int> n[D+1];
    public: Property<PropertyType::ReadOnly, THIS, PropertySetGet<THIS, T>, T> lower[D+1];
    public: Property<PropertyType::ReadOnly, THIS, PropertySetGet<THIS, T>, T> upper[D+1];
    public: Property<PropertyType::ReadOnly, THIS, PropertySet<THIS, T>, T> theta;
    public: HeatEquation() = default;
    public: HeatEquation(T theta, T lower[D+1], T upper[D+1], int n[D+1]) :
        Boundary<T,D>(lower, upper, n),
        theta(PROPERTY(this->theta, theta, Delegate<THIS,T,T>(this, &THIS::Theta))) {

        for(int i = 0; i < D; i++)
            u[i] = new T[0];
        Init<0>(n[0]);
        Init<1>(n[1]);
        Initialize();
        this->theta = theta;
    }
    public: ~HeatEquation() {
        for(int i = 0; i < D; i++) {
            delete matrix[i];
            delete [] u[i];
        }
    }
    private: template<unsigned int DIM> inline ArrayLength<Boundary<T,D>,T> Data() {
        return Boundary<T,D>::n[DIM];
    }
    private: template<unsigned int DIM> void Init(int n) {
        if(DIM) {
            matrix[DIM-1] =  new Matrix<MatrixType::Tridiagonal_m1_C_m1, T>(n-2);
            delete [] u[DIM-1];
            u[DIM-1] = new T[n];
        }
        this->lower[DIM] = Property<PropertyType::ReadOnly, THIS, PropertySetGet<THIS, T>, T>(PropertySetGet<THIS, T>(Delegate<THIS,T,T>(this, &THIS::Lower<DIM>),Delegate<THIS,T,void>(this, &THIS::Lower<DIM>)));
        this->upper[DIM] = Property<PropertyType::ReadOnly, THIS, PropertySetGet<THIS, T>, T>(PropertySetGet<THIS, T>(Delegate<THIS,T,T>(this, &THIS::Upper<DIM>),Delegate<THIS,T,void>(this, &THIS::Upper<DIM>)));
        this->n[DIM] = Property<PropertyType::ReadOnly, THIS, PropertySetGet<THIS, int>, int>(PropertySetGet<THIS, int>(Delegate<THIS,int,int>(this,&THIS::N<DIM>), Delegate<THIS,int,void>(this,&THIS::N<DIM>)));
        this->X[DIM] = Property<PropertyType::ReadOnly, THIS, PropertyGet<THIS, ArrayLength<Boundary<T,D>,T>>, ArrayLength<Boundary<T,D>,T>>(PropertyGet<THIS, ArrayLength<Boundary<T,D>,T>>(Delegate<THIS,ArrayLength<Boundary<T,D>,T>,void>(this,&THIS::Data<DIM>)));
    }
    public: void Initialize() {
        for(int i = 0; i < D; i++) {
            for(int j = 0; j < n[i+1]; j++)
                u[i][j] = 0;
        }
    }
    private: template<unsigned int DIM> inline T Lower() {
        return Boundary<T,D>::lower[DIM];
    }
    private: template<unsigned int DIM> inline T Lower(T lower) {
        Boundary<T,D>::lower[DIM] = lower;
        Theta(theta);
        return lower;
    }
    private: template<unsigned int DIM> inline int N() {
        return  Boundary<T,D>::n[DIM];
    }
    private: template<unsigned int DIM> inline int N(int n) {
        Boundary<T,D>::n[DIM] = n;
        if(DIM) {
            matrix[DIM-1]->n = n-2;
            delete [] u[DIM-1];
            u[DIM-1] = new T[n];
        }
        return n;
    }
    public: void Print(ofstream& file) {
        for(int i = 1; i <= D; i++) {
            if(i > 1)
                file << endl;
            for(int j = 0; j < n[i]; j++) {
                if(j)
                    file << '\t';
                file << this->x[i][j];
            }
        }
        for(int i = 0; i < D; i++) {
            file << endl;
            for(int j = 0; j < n[i+1]; j++) {
                if(j)
                    file << '\t';
                file << this->u[i][j];
            }
        }
    }
    public: void Solve() {
        T* _u[D], *_v[D];           // Preparations
        T* u, *v;
        T alpha;
        int n[D+1], _n;
        n[0] = this->n[0];
        for(int i = 1; i <= D; i++)
            n[i] = this->n[i]-2;
        for(int i = 0; i < D; i++) {
            _v[i] = new T[n[i+1]+2];
            _v[i][0] = this->u[i][0];
            _v[i][n[i+1]+1] = this->u[i][n[i+1]+1];
            _u[i] = this->u[i];
        }

        if(theta) {            // Implicit
            for(int i = 1, j, k; i < n[0]; i++) {
                for(j = 0; j < D; j++) {
                    u = this->u[j];
                    v = _v[j];
                    alpha = this->alpha[j];
                    _n = n[j+1];

                    v[1] = _1_theta*(v[2] - 2*v[1] + v[0]) + alpha*v[1] + v[0];
                    for(k = 2; k < _n; k++)
                        v[k] = _1_theta*(v[k+1] - 2*v[k] - v[k-1]) + alpha*v[k];
                    v[_n] = _1_theta*(v[_n+1] - 2*v[_n] - v[_n-1]) + alpha*v[_n] + v[_n+1];

                    matrix[j]->Solve(&v[1], _n);

                    this->u[j] = v;
                }
            }
        } else {         // Explicit
            for(int i = 1, j, k; i < n[0]; i++) {
                for(j = 0; j < D; j++) {
                    u = this->u[j];
                    v = _v[j];
                    alpha = this->alpha[j];
                    _n = n[j+1];

                    for(k = 1; k <= _n; k++)
                        v[k] = u[k] + alpha*(u[k+1] - 2*u[k] + u[k-1]);

                    this->u[j] = v;
                }
            }

        }

        for(int i = 0; i < D; i++)  // Clean up
            delete [] _u[i];
    }
    private: inline T Theta(T theta) {
        if(theta) {
            _1_theta = 1-(T)1/theta;
            T dt =  (upper[0]-lower[0])/n[0];
            for(int i = 0; i < D; i++) {
                alpha[i] = (upper[i+1]-lower[i+1])/n[i+1];
                alpha[i] *= alpha[i];
                alpha[i] /= dt* theta;
                matrix[i]->Diagonal(2+alpha[i]);
            }
        } else {
            for(int i = 0; i < D; i++) {
                alpha[i] = (upper[i+1]-lower[i+1])/n[i+1];
                alpha[i] *= alpha[i];
                alpha[i] = (upper[0]-lower[0])/(n[0]*alpha[i]);
            }
        }
        return theta;
    }
    private: template<unsigned int DIM> inline T Upper() {
        return Boundary<T,D>::upper[DIM];
    }
    private: template<unsigned int DIM> inline T Upper(T upper) {
        Boundary<T,D>::upper[DIM] = upper;
        Theta(theta);
        return upper;
    }
};

#endif // HEATEQUATION_H
