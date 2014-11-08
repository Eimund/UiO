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
            matrix[DIM-1] =  new Matrix<MatrixType::Tridiagonal_m1_C_m1, T>(n);
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
        if(DIM)
            matrix[DIM-1]->n = n;
        return n;
    }
    private: inline T Theta(T theta) {
        if(theta) {
            _1_theta = 1-(T)1/theta;
            T dt =  this->upper[0]-this->lower[0];
            for(int i = 0; i < D; i++) {
                alpha[i] = (this->upper[i+1]-this->lower[i+1]);
                alpha[i] *= alpha[i];
                alpha[i] /= dt* theta;
                matrix[i]->Diagonal(2+alpha[i]);
            }
        } else {
            for(int i = 0; i < D; i++) {
                alpha[i] = (this->upper[i+1]-this->lower[i+1]);
                alpha[i] *= alpha[i];
                alpha[i] = (this->upper[0]-this->lower[0])/alpha[i];
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
