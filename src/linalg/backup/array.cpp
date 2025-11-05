#include "array.h"

template<class T> 
Array<T>::Array() : _data(0), _size(0)
{}

template<class T> 
Array<T>::Array(uint n) :  _data(new T[n]), _size(n)
{}

template<class T>
Array<T>::Array(const Array<T>& v) :  _data(new T[v._size]), _size(v._size)
{
  for (uint i=0; i<size(); i++)
  {
    (*this)(i) = v(i);
  }
}

/*template<class T> 
Array<T>::~Array()
{

}
*/

template<class T> 
uint Array<T>::size() const
{
  return _size;
}

// with size check
template<class T> 
const T& Array<T>::operator()(uint i) const
{
#ifdef DEBUG
  if (i<0 || i>size())
  {
    std::cerr << "*** error in Array<T>::operator()(uint i) const: index out of range\n";
    exit(0);
  }
#endif
  return _data[i];
}

// with size check
template<class T>
T& Array<T>::operator()(uint i)
{
#ifdef DEBUG
  if (i<0 || i>size())
  {
    std::cerr << "*** error in double Array<T>::operator()(uint i): index out of range\n";
    exit(0);
  }
#endif
  return _data[i];
}

template<class T> 
const T& Array<T>::operator[](uint i) const
{
  return _data[i];
}

template<class T> 
T& Array<T>::operator[](uint i)
{
  return _data[i];
}

template<class T> 
Array<T>& Array<T>::operator=(const Array<T>& v)
{
  resize(v.size());
  for (uint i=0; i<size(); i++)
    (*this)(i) = v(i);

  return *this;
}

template<class T>
Array<T>& Array<T>::operator=(const T& d)
{
  for (uint i=0; i<size(); i++)
    (*this)(i) = d;

  return *this;
}

template<class T>
Array<T>& Array<T>::resize(uint n)
{
  if (_size != n)
  {
    delete [] _data;
    _data = new T[n];
    _size = n;
  }
  return *this;
}

template<class T>
T* Array<T>::ref() const
{
  return _data;
}

//
// I/O
// cout << v << endl;
//
template<class T> 
std::ostream& operator<<(std::ostream& o, const Array<T>& v)
{
  for (uint i=0; i<v.size(); i++)
  {
    o << v(i) << " ";
  }
  return o;
}

template<class T>
std::istream& operator>>(std::istream& is, Array<T>& v)
{
  uint n;
  is >> n;
  v.resize(n);
  for (uint i=0; i<n; i++)
    is >> v(i);

  return is;
}

