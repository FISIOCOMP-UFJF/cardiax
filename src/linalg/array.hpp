#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <iostream>
#include <cstdlib>

using namespace std;

template<class T=double>
class Array
{ 
public:
  Array();
  Array(int n);
  Array(const Array<T>&);
  ~Array();
  
  const T& operator()(int) const;            // element access:   d=v(n)
  T& operator()(int);                        // element access:   v(n)=d
  const T& operator[](int) const;            // element access, no size checks
  T& operator[](int);                        // element access, no size checks
  
  Array<T>& operator=(const Array<T>& v2);   // v=v2;
  Array<T>& operator=(const T& d);           // v=d;  sets all element to d
  Array<T>& resize(int);
  
  T  max();
  T* ref() const { return _data; }
  int size() const { return _size; } 

protected:
  T *_data;
  int _size;
};

// output
template<class T>
std::ostream& operator<<(std::ostream& o, const Array<T>& v);

// input
template<class T>
std::istream& operator>>(std::istream& is, Array<T>& v);

// Implementation -------------------------------------------------------------

template<class T> Array<T>::Array()
  : _data(0), _size(0)
{}

template<class T> Array<T>::Array(int n)
  :  _data(new T[n]), _size(n)
{}

template<class T> Array<T>::Array(const Array<T>& v)
  :  _data(new T[v._size]), _size(v._size)
{
  for (int i=0; i<size(); i++)
  {
    (*this)(i) = v(i);
  }
}

template<class T> Array<T>::~Array()
{
  delete [] _data;
}

template<class T> T Array<T>::max()
{
  T m;
  m = _data[0];
  for(int i=1; i<_size; i++)
    if (_data[i] > m) m = _data[i];
  return m;
}

// ----------------------------------------------------------------------------

template<class T> const T& Array<T>::operator()(int i) const
{
#ifdef DEBUG
  if (i<0 || i>size())
    throw std::runtime_error("Array<T>::operator()(int): index out of range");
#endif
  return _data[i];
}

template<class T> T& Array<T>::operator()(int i)
{
#ifdef DEBUG
  if (i<0 || i>size())
    throw std::runtime_error("Array<T>::operator()(int): index out of range");
#endif
  return _data[i];
}

template<class T> const T& Array<T>::operator[](int i) const
{
  return _data[i];
}

template<class T> T& Array<T>::operator[](int i)
{
  return _data[i];
}

// ----------------------------------------------------------------------------

template<class T> Array<T>& Array<T>::operator=(const Array<T>& v)
{
  resize(v.size());
  for (int i=0; i<size(); i++)
  {
    (*this)(i) = v(i);
  }
  return *this;
}

template<class T> Array<T>& Array<T>::operator=(const T& d)
{
  for (int i=0; i<size(); i++)
  {
    (*this)(i) = d;
  }
  return *this;
}


template<class T> Array<T>& Array<T>::resize(int n)
{
  if (_size != n)
  {
    delete [] _data;
    _data = new T[n];
    _size = n;
  }
  return *this;
}

// output: cout << v << endl;
template<class T> std::ostream& operator<<(std::ostream& o, const Array<T>& v)
{
  for (int i=0; i<v.size(); i++)
  {
    o << v(i) << " ";
  }
  return o;
}

template<class T> std::istream& operator>>(std::istream& is, Array<T>& v)
{
  int n;
  is >> n;
  v.resize(n);
  for (int i=0; i<n; i++)
    is >> v(i);

  return is;
}

#endif
