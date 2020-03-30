/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Author: Sanjib Sharma                                                     *
 * Copyright (c) 2012 Sanjib Sharma                                          *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of Galaxia. The full Galaxia copyright notice, including*
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and COPYRIGHT, which can be found at the root           *
 * of the source code distribution tree.                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef MATRIX_H_
#define MATRIX_H_

#include "sutils.h"

template <class T>
class Matrix_base {
protected:
	int nn;
	int mm;
	T **v;
public:
	Matrix_base(): nn(0), mm(0), v(NULL) {}
	// Zero-based array
	Matrix_base(int n, int m): nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
	{
		int i, nel = m * n;
		if (v)
			v[0] = nel > 0 ? new T[nel] : NULL;
		for (i = 1; i < n; i++)
			v[i] = v[i - 1] + m;
	}
	//Initialize to constant
	Matrix_base(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
	{
		int i, j, nel = m * n;
		if (v)
			v[0] = nel > 0 ? new T[nel] : NULL;
		for (i = 1; i < n; i++)
			v[i] = v[i - 1] + m;
		for (i = 0; i < n; i++)
			for (j = 0; j < m; j++)
				v[i][j] = a;
	}
	// Initialize to array
	Matrix_base(int n, int m, const T *a)	: nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
	{
		int i,j,nel=m*n;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< n; i++) v[i] = v[i-1] + m;
		for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
	}
	// Copy constructor
	Matrix_base(const Matrix_base &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn > 0 ? new T*[nn] : NULL)
	{
		int i, j, nel = mm * nn;
		if (v)
			v[0] = nel > 0 ? new T[nel] : NULL;
		for (i = 1; i < nn; i++)
			v[i] = v[i - 1] + mm;
		for (i = 0; i < nn; i++)
			for (j = 0; j < mm; j++)
				v[i][j] = rhs[i][j];
	}
	//assignment
	// postcondition: normal assignment via copying has been performed;
	//		if matrix and rhs were different sizes, matrix
	//		has been resized to match the size of rhs
	Matrix_base & operator=(const Matrix_base &rhs)
	{
		if (this != &rhs)
		{
			int i,j,nel;
			if (nn != rhs.nn || mm != rhs.mm) {
				if (v != NULL) {
					delete[] (v[0]);
					delete[] (v);
				}
				nn=rhs.nn;
				mm=rhs.mm;
				v = nn>0 ? new T*[nn] : NULL;
				nel = mm*nn;

				if (v) v[0] = nel>0 ? new T[nel] : NULL;
				for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
			}
			for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
		}
		return *this;
	}
	// make T available externally
	typedef T value_type;
	inline T* operator[](const int i){	return v[i];}	//subscripting: pointer to row i
	inline const T* operator[](const int i) const{	return v[i];}
	inline int nrows() const {return nn;}
	inline int ncols() const {return mm;}
	~Matrix_base()
	{
		if (v != NULL)
		{
			delete[] (v[0]);
			delete[] (v);
		}
	}
};














template <class T>
class Matrix {
private:
	int nn;
	int mm;
	T **v;
public:
	Matrix();
	Matrix(int n, int m);			// Zero-based array
	Matrix(int n, int m, const T &a);	//Initialize to constant
	Matrix(int n, int m, const T *a);	// Initialize to array
	Matrix(const Matrix &rhs);		// Copy constructor
	Matrix & operator=(const Matrix &rhs);	//assignment
	typedef T value_type; // make T available externally
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline void identity()
	{
		for(int j=0;j<nn;++j)
			for(int k=0;k<mm;++k)
			{
				if(k==j)
					v[j][k]=1.0;
				else
					v[j][k]=0.0;
			}

	}
	inline void transpose()
	{
		T temp;
		assert(nn=mm);
		for(int j=0;j<nn;++j)
			for(int k=j+1;k<mm;++k)
			{
				temp=v[k][j]; v[k][j]=v[j][k];v[j][k]=temp;
			}

	}
	Matrix& operator*=(const T &a)
	{
		for(int j=0;j<nn;++j)
			for(int k=0;k<mm;++k)
				v[j][k]*=a;
		return *this;

	}
	Matrix& operator*=(const Matrix &M1)	//multiply
	{
		assert(nn=mm);
		Matrix_base<T> temp(nn,mm);
		for(int j=0;j<nn;++j)
			for(int k=0;k<mm;++k)
			{
				temp[j][k]=v[j][k];
				v[j][k]=0.0;
			}
		for(int j=0;j<nn;++j)
			for(int k=0;k<mm;++k)
				for(int l=0;l<mm;++l)
					v[j][k]+=temp[j][l]*M1[l][k];
		return *this;
	}
	inline int nrows() const;
	inline int ncols() const;
	void resize(int newn, int newm); // resize (contents not preserved)
	void assign(int newn, int newm, const T &a); // resize and assign a constant value
	void print()
	{
		std::cout<<"Size: "<<nn<<"x"<<mm<<std::endl;
		for(int j=0;j<nn;++j)
		{
			for(int k=0;k<mm;++k)
			{
				std::cout<<v[j][k]<<" ";
			}
			std::cout<<std::endl;
		}
	}
	template <class T1,class T2>
	void mult(T1 x,T2 y)
	{
		vector<T> temp(nn);
		for(int j=0;j<nn;++j)
		{
			temp[j]=x[j];
			y[j]=0.0;
		}
		for(int j=0;j<nn;++j)
			for(int k=0;k<mm;++k)
				y[j]+=v[j][k]*temp[k];
	}
	template <class T1,class T2>
	void Tmult(T1 x,T2 y)
	{
		vector<T> temp(nn);
		for(int j=0;j<nn;++j)
		{
			temp[j]=x[j];
			y[j]=0.0;
		}
		for(int j=0;j<nn;++j)
			for(int k=0;k<mm;++k)
				y[j]+=v[k][j]*temp[k];
	}
	~Matrix();
};

template <class T>
Matrix<T>::Matrix() : nn(0), mm(0), v(NULL) {}

template <class T>
Matrix<T>::Matrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
Matrix<T>::Matrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
Matrix<T>::Matrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
Matrix<T>::Matrix(const Matrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
	int i,j,nel=mm*nn;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
Matrix<T> & Matrix<T>::operator=(const Matrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		int i,j,nel;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = nn>0 ? new T*[nn] : NULL;
			nel = mm*nn;

			if (v) v[0] = nel>0 ? new T[nel] : NULL;
			for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
		}
		for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
inline T* Matrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* Matrix<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline int Matrix<T>::nrows() const
{
	return nn;
}

template <class T>
inline int Matrix<T>::ncols() const
{
	return mm;
}

template <class T>
void Matrix<T>::resize(int newn, int newm)
{
	int i,nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
}

template <class T>
void Matrix<T>::assign(int newn, int newm, const T& a)
{
	int i,j,nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
}

template <class T>
Matrix<T>::~Matrix()
{
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}




#endif /* MATRIX_H_ */
