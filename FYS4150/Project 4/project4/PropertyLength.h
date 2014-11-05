#ifndef PROPERTYLENGTH_H
#define PROPERTYLENGTH_H

template<class T> class PropertyLength {
    private: T** arr;
    private: unsigned int length;
    public: PropertyLength(T* &arr, unsigned int length) : arr(arr) {
        length = 1;
        *arr = new T[1];
    }
    public: ~PropertyLength() {
        delete [] *arr;
    }
    public: unsigned int& operator = (const unsigned int& length) {
        if(length && length != this->length) {
            T* arr = new T[length];
            unsigned int n = length < this->length ? length : this->length;
            for(unsigned int i = 0; i < n; i++)
                arr[i] = (*this->arr)[i];
            delete [] *this->arr;
            *this->arr = arr;
            return this->length = length;
        } else
            return this->length;
    }
    public: operator unsigned int () const {
        return length;
    }
};

#endif // PROPERTYLENGTH_H
