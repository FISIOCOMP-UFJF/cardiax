#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <cstdlib>
#include "datatype.h"

template<class T> 
class Array
{
  private:
    uint _size;
    T *_data;

  public:
    Array();
    Array(uint n);
    Array(const Array<T>&);
    ~Array(){
      delete [] _data;
    }

    uint size() const;

    const T& operator()(uint) const;  // element access:   d=v(n)
    T& operator()(uint);              // element access:   v(n)=d
    const T& operator[](uint) const;  // element access, no size checks
    T& operator[](uint);              // element access, no size checks

    Array<T>& operator=(const Array<T>& v2); // v=v2;
    Array<T>& operator=(const T& d);         // v=d;  sets all element to d
    Array<T>& resize(uint);

    T* ref() const;                          // return reference to storage
};

#endif
