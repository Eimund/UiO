#ifndef GRAVITATION_H
#define GRAVITATION_H

#include <cmath>
#include <fstream>
#include <string>
#include "Differential.h"

using namespace std;

enum class CircularDirection {
    CW,
    CCW
};

template<class T> struct Null {
    static const T value;
};
template<class T> const T Null<T>::value = 0;
template<class T> struct Null<T*> {
    static const T* value;
};
template<class T> const T* Null<T*>::value = nullptr;

template<class T> struct Array {
    T element;
    Array<T>* prev;
    Array<T>* next;
    Array() {
        prev = this;
        next = this;
    }
    ~Array() {
        if(next != this)
            delete next;
    }
    T operator[] (const int i) {
        if(prev != this) {
            if(i) {
                if(i > 0)
                    return (*next)[i-1];
                else
                    return (*prev)[i+1];
            }
            return this->element;
        } else if(next != this)
            return (*next)[i];
        return (T)Null<T>::value;
    }
    Array<T>* Index (const int i) {
        if(prev != this) {
            if(i) {
                if(i > 0)
                    return next->Index(i-1);
                else
                    return prev->Index(i+1);
            }
            return this;
        } else if(next != this)
            return next->Index(i);
        return nullptr;
    }
    void Add(T element) {
        if(next != this)
            next->Add(element);
        else {
            next = new Array<T>;
            next->prev = this;
            next->next = next;
            next->element = element;
        }
    }
    int Length() {
        if(next != this)
            return next->Length()+1;
        return 0;
    }
};
template<class T, unsigned int DIM> struct Body {
    string name;
    T m;
    T x0[DIM], v0[DIM];
    T* x[DIM];
    T* v[DIM];
    Body(string name, T m, T x0[DIM], T v0[DIM], unsigned int len) {
        this->name = name;
        this->m = m;
        for(unsigned int i = 0; i < DIM; i++) {
            x[i] = new T[len];
            v[i] = new T[len];
            x[i][0] = x0[i];
            v[i][0] = v0[i];
            this->x0[i] = x0[i];
            this->v0[i] = v0[i];
        }
    }
    ~Body() {
        for(unsigned int i = 0; i < DIM; i++) {
            delete [] x[i];
            delete [] v[i];
        }
    }
    void Print(ofstream& file, unsigned int len) {
        file << name;
        for(unsigned int j = 0, i; j < DIM; j++) {
            file << endl;
            for(i = 0; i < len; i++) {
                if(i)
                    file << '\t';
                file << x[j][i];
            }
        }
        for(unsigned int j = 0, i; j < DIM; j++) {
            file << endl;
            for(i = 0; i < len; i++) {
                if(i)
                    file << '\t';
                file << v[j][i];
            }
        }
    }
    void SetLength(unsigned int len) {
        for(unsigned int i = 0; i < DIM; i++) {
            delete [] x[i];
            delete [] v[i];
            x[i] = new T[len];
            v[i] = new T[len];
            x[i][0] = x0[i];
            v[i][0] = v0[i];
        }
    }
};
template<class T, unsigned int DIM> class System : Differential_2<T> {
    protected: unsigned int _length;
    private: class _length_ {
        private: System<T, DIM>* owner;
        public: _length_(System<T, DIM>* owner) : owner(owner) {
            this->owner->_length = 1;
        }
        public: unsigned int& operator = (const unsigned int& length) {
            if(owner->_length != length) {
                auto body = owner->body;
                unsigned int j = 0;
                delete [] owner->t;
                owner->t = new T[length];
                while(body->next != body) {
                    body = body->next;
                    body->element->SetLength(length);
                    for(unsigned int i = 0; i < DIM; i++, j++) {
                        owner->x[j] = body->element->x[i];
                        owner->v[j] = body->element->v[i];
                    }
                }
                return owner->_length = length;
            }
            return owner->_length;
        }
        public: operator unsigned int () const {
            return owner->_length;
        }
    };
    public: _length_ length;
    private: unsigned int width;
    public: T G;
    private: T* t, *a;
    private: T** x, **v;
    private: Array<Body<T,DIM>*>* body;
    public: System(T G) : length(this) {
        width = 0;
        this->G = G;
        t = new T[_length];
        a = new T[0];
        x = new T*[0];
        v = new T*[0];
        body = new Array<Body<T,DIM>*>;
    }
    public: ~System() {
        delete [] a;
        delete [] x;
        delete [] v;
        delete [] t;
        delete body;
    }
    public: void Add(string name, T m, T x0[DIM], T v0[DIM]) {
        body->Add(new Body<T,DIM>(name, m, x0, v0, _length));

        unsigned int len = body->Length();
        width = DIM*len;
        delete [] a;
        delete [] x;
        delete [] v;
        a = new T[width];
        x = new T*[width];
        v = new T*[width];

        for(unsigned int i = 0, j, k = 0; i < len; i++, k+=DIM) {
            for(j = 0; j < DIM; j++) {
                x[k+j] = (*body)[i]->x[j];
                v[k+j] = (*body)[i]->v[j];
            }
        }
    }
    public: template<unsigned int x, unsigned int y> T* Coplanar(Body<T,DIM>* , Body<T,DIM>* ) {
        T* alpha = new T[DIM-2];
        return alpha;
    }
    public: T Distance(unsigned int i, Body<T,DIM>* b1, Body<T,DIM>* b2) {
        T diff, dis = 0;
        for(unsigned int j = 0; j < DIM; j++) {
            diff = b1->x[j][i]-b2->x[j][i];
            dis += diff*diff;
        }
        return sqrt(dis);
    }
    public: Body<T,DIM>* Find(string body) {
        auto b = this->body;
        do {
            b = b->next;
            if(!b->element->name.compare(body))
                return b->element;

        } while(b != b->next);
        return nullptr;
    }
    public: template<unsigned int x, unsigned int y, CircularDirection D> void GoesAround(string satellite, string primary) {
        auto sat = Find(satellite);
        auto pri = Find(primary);
        T* alpha = Coplanar<x,y>(sat,pri);
        T r = Distance(0, sat, pri);
        T v = sqrt(G*pri->m/r);
        if(D == CircularDirection::CW) {
            sat->v[x][0] = v*(sat->x[y][0]-pri->x[y][0])/r + pri->v[x][0];
            sat->v[y][0] = v*(pri->x[x][0]-sat->x[x][0])/r + pri->v[y][0];
        } else {
            sat->v[x][0] = v*(pri->x[y][0]-sat->x[y][0])/r + pri->v[x][0];
            sat->v[y][0] = v*(sat->x[x][0]-pri->x[x][0])/r + pri->v[x][0];
        }
        sat->v0[x] = sat->v[x][0];
        sat->v0[y] = sat->v[y][0];
        for(unsigned int i = 0; i < DIM-2; i++)
            alpha[i] *= -1;
        Rotate<x,y>(sat,pri,alpha);
        delete [] alpha;
    }
    private: T* Gravity(T* x) {
        T d,r, diff[DIM];
        auto b1 = body->next;
        decltype(b1) b2;
        for(unsigned int i = 0; i < width; i++)
            a[i] = 0;
        for(unsigned int i = 0, j, k, i1, j1; i < width; i+=DIM, b1=b1->next) {
            for(j=i+DIM, b2=b1->next; j < width; j+=DIM, b2=b2->next) {
                r = 0;
                for(k=0, i1=i, j1=j; k < DIM; k++, i1++, j1++) {
                    diff[k] = x[i1]-x[j1];
                    r += diff[k]*diff[k];
                }
                if(r) {
                    r = G/(r*sqrt(r));
                    for(k=0, i1=i, j1=j; k < DIM; k++, i1++, j1++) {
                        d = r*diff[k];
                        a[i1] -= d*b2->element->m;
                        a[j1] += d*b1->element->m;
                    }
                }
            }
        }
        return a;
    }
    public: void Print(ofstream& file) {
        for(unsigned int i = 0; i < _length; i++) {
            if(i)
                file << '\t';
            file << t[i];
        }
        file << endl;

        auto body = this->body;
        while(body != body->next) {
            if(body != this->body)
                file << endl;
            body = body->next;
            body->element->Print(file, _length);
        }
    }
    public: void PrintEnergyMomentum(ofstream& file) {
        for(unsigned int i = 0; i < _length; i++) {
            if(i)
                file << '\t';
            file << t[i];
        }
        file << endl;

        T value;
        for(unsigned int i = 0; i < _length; i++) {
            value = 0;
            for(unsigned int j = 0; j < body->Length(); j++) {
                for(unsigned int k = 0; k < DIM; k++)
                    value += 0.5*(*body)[j]->m*(*body)[j]->v[k][i]*(*body)[j]->v[k][i];
            }
            if(i)
                file << '\t';
            file << value;
        }
        file << endl;

        for(unsigned int i = 0; i < _length; i++) {
            for(unsigned int j = 0; j < body->Length()-1; j++) {
                value = 0;
                for(unsigned int k = j+1; k < body->Length(); k++)
                    value -= G*(*body)[j]->m*(*body)[k]->m/Distance(i,(*body)[j], (*body)[k]);
                if(i)
                    file << '\t';
                file << value;
            }
        }


        for(unsigned int k = 0; k < DIM; k++) {
            file << endl;
            for(unsigned int i = 0; i < _length; i++) {
                value = 0;
                for(unsigned int j = 0; j < body->Length(); j++)
                    value += (*body)[j]->m*(*body)[j]->v[k][i];
                if(i)
                    file << '\t';
                file << value;
            }
        }
    }
    public: template<unsigned int x, unsigned int y> void Rotate(Body<T,DIM>* , Body<T,DIM>* , T* ) {
    }
    public: template<DifferentialType Type> void Run(T t, unsigned int n) {
        T dt = t/(n-1);
        length = n;
        this->template Solve<Type>(this, &System<T,DIM>::Gravity, this->t, x, v, dt, width, n);
    }
    public: void ZeroMomentum(string body) {
        auto b = Find(body);
        for(unsigned int i = 0; i < DIM; i++) {
            b->v0[i] = 0;
            auto s = this->body;
            do {
                s = s->next;
                if(s->element != b)
                    b->v0[i] -= s->element->m*s->element->v0[i];

            } while(s != s->next);
            b->v0[i] /= b->m;
            b->v[i][0] = b->v0[i];
        }
    }
};

#endif // GRAVITATION_H
