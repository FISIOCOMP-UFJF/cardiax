#ifndef _MALLOC_H_
#define _MALLOC_H_

/**
 * @file   malloc.h
 * @author Bernardo M. Rocha <bernardo@lncc.br>
 * @date   Wed Mar  3 11:05:21 2010
 * 
 * @brief  Dynamic memory allocation header file
 * 
 * 
 */

#include <cstdio>
#include <cstdlib>

using namespace std;

template <typename T>
T * new_array(const size_t N)
{
  return (T*) calloc(N, sizeof(T));
}

template <typename T>
T * resize_array(T * v, size_t new_size)
{
  v = (T *) realloc(v, new_size*sizeof(T));
  if(v == NULL){
    cerr << "resize_array: error reallocating memory!" << endl;
    exit(1);
  }
  return v;
}

template <typename T>
void delete_array(T * ptr)
{
  free(ptr);
}

template <typename T>
T ** new_2D_array(const size_t M, const size_t N)
{
  size_t i;
  T **mat;
  mat = (T **) calloc(M, sizeof(T*));
  if(mat == NULL) fprintf(stderr,"allocation error 1 in dmatrix");
  for(i=0; i<M; i++){
    mat[i] = (T *) calloc(N, sizeof(T));
    if(mat[i] == NULL) fprintf(stderr,"allocation error 2 in dmatrix");
  }
  return mat;
}

template <typename T>
void delete_2D_array(T ** ptr, const size_t M, const size_t N)
{
  size_t i;
  for(i=0; i<M; i++)
    free(ptr[i]);

  free(ptr);
}

template <typename T>
T *** new_3D_array(const size_t M, const size_t N, const size_t P)
{
  size_t i,j;
  T ***mat;
  mat = (T ***) calloc(M, sizeof(T**));
  if(mat == NULL) 
    fprintf(stderr,"allocation error 1 in dmatrix");
  for(i=0; i<M; i++)
  {
    mat[i] = (T **) calloc(N, sizeof(T*));
    if(mat[i] == NULL) 
      fprintf(stderr,"allocation error 2 in dmatrix");
    for(j=0; j<N; j++)
    {
      mat[i][j] = (T *) calloc(P, sizeof(T));
      if(mat[i][j] == NULL) 
	fprintf(stderr,"allocation error 3 in dmatrix");
    }
  }
  return mat;
}

template <typename T>
void delete_3D_array(T *** ptr, const size_t M, 
		     const size_t N, const size_t P)
{
  size_t i,j;
  for(i=0; i<M; i++){
    for(j=0; j<N; j++){
      free(ptr[i][j]);
    }
    free(ptr[i]);
  }
  free(ptr);
}

#endif
