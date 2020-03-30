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

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_
#include <algorithm>
#include<iterator>
#include "Interp.h"

//using namespace std;

struct Ran {
	uint64_t u,v,w;
	Ran(uint64_t j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	inline uint64_t int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		uint64_t x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	inline uint32_t int32() { return (uint32_t)int64(); }
};
struct Normaldev : Ran {
	double mu,sig;
	Normaldev(double mmu, double ssig, uint64_t i)
	: Ran(i), mu(mmu), sig(ssig){}
	double dev() {
		double u,v,x,y,q;
		do {
			u = doub();
			v = 1.7156*(doub()-0.5);
			x = u - 0.449871;
			y = fabs(v) + 0.386595;
			q = x*x + y*(0.19600*y-0.25472*x);
		} while (q > 0.27597
			&& (q > 0.27846 || (v*v) > -4.*log(u)*(u*u)));
		return mu + sig*v/u;
	}
};

//struct Ran
//{
//	unsigned long long u, v, w;
//	Ran(unsigned long long j) :		v(4101842887655102017LL), w(1)
//	{
//		u = j ^ v;
//		int64();
//		v = u;
//		int64();
//		w = v;
//		int64();
//	}
//
//	inline unsigned long long int64()
//	{
//		u = u * 2862933555777941757LL + 7046029254386353087LL;
//		v ^= v >> 17;
//		v ^= v << 31;
//		v ^= v >> 8;
//		w = 4294957665U * (w & 0xffffffff) + (w >> 32);
//		unsigned long long x = u ^ (u << 21);
//		x ^= x >> 35;
//		x ^= x << 4;
//		return (x + v) ^ w;
//	}
//
//	inline double doub()
//	{
//		return 5.42101086242752217E-20 * int64();
//	}
//
//	inline unsigned int int32()
//	{
//		return (unsigned int) int64();
//	}
//};

//class Normaldev: Ran
//{
//public:
//	double mu, sig;
//	Normaldev(double mmu, double ssig, unsigned long long i):Ran(i), mu(mmu), sig(ssig)
//	{
//	}
//	double dev()
//	{
//		double u, v, x, y, q;
//		do
//		{
//			u = doub();
//			v = 1.7156 * (doub() - 0.5);
//			x = u - 0.449871;
//			y = abs(v) + 0.386595;
//			q = (x*x) + y * (0.19600 * y - 0.25472 * x);
//		} while (q > 0.27597 && (q > 0.27846 || (v*v) > -4. * log(u) * (u*u)));
//		return mu + sig * v / u;
//	}
//};


extern Ran nrRan;
extern Normaldev nrGauss;

//static Ran nrRan(13);
//static Normaldev nrGauss(0.0,1.0,17);
//static Ran nrRan(33);
//static Normaldev nrGauss(0.0,1.0,33);

inline double grandomu()
{
//	return drand48();
	return nrRan.doub();
}

inline double grandomn()
{
	return nrGauss.dev();
}



template<class T>
double theta(T Pos)
{
	return acos(sqrt((Pos[0]*Pos[0]+Pos[1]*Pos[1])/(Pos[0]*Pos[0]+Pos[1]*Pos[1]+Pos[2]*Pos[2])));
}

template<class T,class T1>
double distance(T Pos,T1 PosC,int j)
{
	double dist=0.0;
	for(int i=0;i<j;++i)
		dist+=(Pos[i]-PosC[i])*(Pos[i]-PosC[i]);
	return sqrt(dist);
}

template<class T>
void printv(T a,T b)
{
	std::cout<<"Printing Vector"<<endl;
	while(a!=b)
//		std::cout<<*a++<<" "<<endl;
		std::cout<<scientific<<*a++<<endl;
//		std::cout<<scientific<<(*a)<<setw(12)<<endl;
	std::cout<<endl;
}

template<class T>
void Stats(T a,T b)
{
	int ntot=b-a;
	double xmin,xmax,xmean,xstddev;
	std::cout<<"Printing Stats min max mean stddev"<<endl;
	xmin=xmax=xmean=xstddev=(*a);
	a++;
	while(a!=b)
	{
		xmean+=(*a);
		xstddev+=(*a)*(*a);
		if((*a)<xmin)
			xmin=(*a);
		else if((*a)>xmax)
			xmax=(*a);
		a++;
//		std::cout<<*a++<<" "<<endl;
//		std::cout<<scientific<<(*a)<<setw(12)<<endl;
	}
	xmean/=ntot;
	xstddev=sqrt(xstddev/ntot-xmean*xmean);
	std::cout<<scientific<<ntot<<" "<<xmin<<" "<<xmax<<" "<<xmean<<" "<<xstddev<<std::endl;
}


template<class T>
void printv(T a,T b,const char *s)
{
	std::cout<<"Printing Vector "<<s<<endl;
	while(a!=b)
		std::cout<<*a++<<" "<<endl;
//		std::cout<<scientific<<*a<<endl;
	std::cout<<endl;
}

template<class T>
void printv(vector<T> &a,const char *s)
{
	std::cout<<"Printing Vector "<<s<<endl;
    for(unsigned int i=0;i<a.size();++i)
		std::cout<<a[i]<<" "<<endl;
//		std::cout<<scientific<<a[i]<<endl;
	std::cout<<endl;
}


template<class T>
double int_tabulated(vector<T> &x,vector<T> &y,T xmin,T xmax)
{
    double temp=0;
    for(unsigned int i=1;i<x.size();++i)
    {
    	if(((x[i])<=xmax)&&(x[i-1]>=xmin))
    		temp+=(x[i]-x[i-1])*(y[i]+y[i-1])*0.5;
    }
    return temp;
}

template<class T>
double int_tabulated(vector<T> &x,vector<T> &y)
{
	double temp=0;
	unsigned int i;
    for(i=1;i<x.size();++i)
    {
	    temp+=(x[i]-x[i-1])*(y[i]+y[i-1])*0.5;
    }
    return temp;
}

template<class T>
void linspace(vector<T>&x,double xmin,double xmax,double h)
{
	int size1=int((xmax-xmin)/h+0.5)+1;
	if(size1 >1)
	{
		x.resize(size1);
		double temp=(xmax-xmin)/(x.size()-1);
		x.front()=xmin;
		for(unsigned int i=1;i<x.size();++i)
			x[i]=xmin+i*temp;
	}
	else
	{
		size1=2;
		x.resize(size1);
		x[0]=xmin;
		x[1]=xmax;
	}
}


template<class T>
void linspace(vector<T>&x,double xmin,double xmax,int size1)
{
	if(size1<2)
	{
		cout<<"linspace: size should be greater than 1"<<endl;
		exit(1);
	}
	x.resize(size1);
	double temp=(xmax-xmin)/(x.size()-1);
	for(unsigned int i=0;i<x.size();++i)
		x[i]=xmin+i*temp;

}


template<class T>
void logspace(vector<T> &x,double xmin,double xmax,int size1)
{
	if(size1<2)
	{
		cout<<"logspace: size should be greater than 1"<<endl;
		exit(1);
	}
	x.resize(size1);
	double temp=(log10(xmax)-log10(xmin))/(x.size()-1);
	x.front()=xmin;
	for(unsigned int i=1;i<x.size();++i)
		x[i]=pow(10.0,log10(xmin)+i*temp);
}

template<class T>
void logspace(vector<T>&x,double xmin,double xmax,double h)
{
	int size1=int((log10(xmax)-log10(xmin))/h+0.5)+1;
	if(size1 >1)
	{
		x.resize(size1);
		double temp=(log10(xmax)-log(xmin))/(x.size()-1);
		x.front()=xmin;
		for(unsigned int i=1;i<x.size();++i)
			x[i]=pow(10.0,log10(xmin)+i*temp);
	}
	else
	{
		size1=2;
		x.resize(size1);
		x[0]=xmin;
		x[1]=xmax;
	}
}

template<class T>
void loglinspace(vector<T>&x,double xmin,double xmax,int size1)
{

	if(size1 >1)
	{
		T xmax1,xmin1;
		int size2=ceil(log10(xmax/xmin))/log10(40.0);
		if(size2<1)
			size2=1;
		x.resize(size2*size1);
		xmax1=xmin*40;
		xmin1=xmin;
		for(int j=0;j<size2;++j)
		{
			if(xmax1>xmax) xmax1=xmax;
			double temp=(xmax1-xmin1)/(size1-1);
			for(int i=0;i<size1;++i)
				x[j*size1+i]=xmin1+i*temp;
			xmin1=xmax1;
			xmax1*=20;
		}
		x.front()=xmin;
		x.back()=xmax;
	}
	else
	{
		size1=2;
		x.resize(size1);
		x[0]=xmin;
		x[1]=xmax;
	}
}



template<class T>
void sort2(T arr, T arr1, T brr)
{
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50
	uint64_t n=arr1-arr;
	arr=arr-1;
	brr=brr-1;

	uint64_t i,ir=n,j,k,l=1;
        int *istack,jstack=0;
        double a,b,temp;

        istack=new int[1+NSTACK];
        for (;;) {
                if (ir-l < M) {
                        for (j=l+1;j<=ir;j++) {
                                a=arr[j];
                                b=brr[j];
                                for (i=j-1;i>=l;i--) {
                                        if (arr[i] <= a) break;
                                        arr[i+1]=arr[i];
                                        brr[i+1]=brr[i];
                                }
                                arr[i+1]=a;
                                brr[i+1]=b;
                        }
                        if (!jstack) {
                                delete[] istack;
                                return;
                        }
                        ir=istack[jstack];
                        l=istack[jstack-1];
                        jstack -= 2;
                } else {
                        k=(l+ir) >> 1;
                        SWAP(arr[k],arr[l+1])
                        SWAP(brr[k],brr[l+1])
                        if (arr[l] > arr[ir]) {
                                SWAP(arr[l],arr[ir])
                                SWAP(brr[l],brr[ir])
                        }
                        if (arr[l+1] > arr[ir]) {
                                SWAP(arr[l+1],arr[ir])
                                SWAP(brr[l+1],brr[ir])
                        }
                        if (arr[l] > arr[l+1]) {
                                SWAP(arr[l],arr[l+1])
                                SWAP(brr[l],brr[l+1])
                        }
                        i=l+1;
                        j=ir;
                        a=arr[l+1];
                        b=brr[l+1];
                        for (;;) {
                                do i++; while (arr[i] < a);
                                do j--; while (arr[j] > a);
                                if (j < i) break;
                                SWAP(arr[i],arr[j])
                                SWAP(brr[i],brr[j])
                        }
                        arr[l+1]=arr[j];
                        arr[j]=a;
                        brr[l+1]=brr[j];
                        brr[j]=b;
                        jstack += 2;
                        if (jstack > NSTACK) exit(1);
                        if (ir-i+1 >= j-l) {
                                istack[jstack]=ir;
                                istack[jstack-1]=i;
                                ir=j-1;
                         } else {
                                 istack[jstack]=j-1;
                                 istack[jstack-1]=l;
                                 l=i;
                         }
                 }
         }
#undef M
#undef NSTACK
#undef SWAP
 }


template<class T,class T1>
void sorti2(T arr, T arr1, T1 brr)
{
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define SWAPI(a1,b1) tempi=(a1);(a1)=(b1);(b1)=tempi;
#define M 7
#define NSTACK 50
	uint64_t n=arr1-arr;
	arr=arr-1;
	brr=brr-1;
	uint64_t i,ir=n,j,k,l=1;
        int *istack,jstack=0;
        double a,temp;
        int b,tempi;

        istack=new int[1+NSTACK];
        for (;;) {
                if (ir-l < M) {
                        for (j=l+1;j<=ir;j++) {
                                a=arr[j];
                                b=brr[j];
                                for (i=j-1;i>=l;i--) {
                                        if (arr[i] <= a) break;
                                        arr[i+1]=arr[i];
                                        brr[i+1]=brr[i];
                                }
                                arr[i+1]=a;
                                brr[i+1]=b;
                        }
                        if (!jstack) {
                                delete[] istack;
                                return;
                        }
                        ir=istack[jstack];
                        l=istack[jstack-1];
                        jstack -= 2;
                } else {
                        k=(l+ir) >> 1;
                        SWAP(arr[k],arr[l+1])
                        SWAPI(brr[k],brr[l+1])
                        if (arr[l] > arr[ir]) {
                                SWAP(arr[l],arr[ir])
                                SWAPI(brr[l],brr[ir])
                        }
                        if (arr[l+1] > arr[ir]) {
                                SWAP(arr[l+1],arr[ir])
                                SWAPI(brr[l+1],brr[ir])
                        }
                        if (arr[l] > arr[l+1]) {
                                SWAP(arr[l],arr[l+1])
                                SWAPI(brr[l],brr[l+1])
                        }
                        i=l+1;
                        j=ir;
                        a=arr[l+1];
                        b=brr[l+1];
                        for (;;) {
                                do i++; while (arr[i] < a);
                                do j--; while (arr[j] > a);
                                if (j < i) break;
                                SWAP(arr[i],arr[j])
                                SWAPI(brr[i],brr[j])
                        }
                        arr[l+1]=arr[j];
                        arr[j]=a;
                        brr[l+1]=brr[j];
                        brr[j]=b;
                        jstack += 2;
                        if (jstack > NSTACK) exit(1);
                        if (ir-i+1 >= j-l) {
                                istack[jstack]=ir;
                                istack[jstack-1]=i;
                                ir=j-1;
                         } else {
                                 istack[jstack]=j-1;
                                 istack[jstack-1]=l;
                                 l=i;
                         }
                 }
         }
#undef M
#undef NSTACK
#undef SWAP
 }




template <class ForwardIterator,class ForwardIterator1>
ForwardIterator unique2 ( ForwardIterator first, ForwardIterator last,ForwardIterator1 first1 )
{
  ForwardIterator result=first;
  ForwardIterator1 result1=first1;
  ++first1;
  while (++first != last)
  {
    if (!(*result == *first))  // or: if (!pred(*result,*first)) for the pred version
    {
    	*(++result)=*first;
    	*(++result1)=*first1;
    }
    first1++;
  }
  return ++result;
}

template <class InputIterator,class T>
T total(InputIterator first, InputIterator last,T val)
{
  while (first != last)
	  val+=*first++;
  return val;
}

template <class T>
T total(vector<T> &v)
{
	T val=0;
	for(unsigned int i=0;i<v.size();++i)
		val+=v[i];
	return val;
}

template <class T1,class T2>
void totalCumul(vector<T1> &v1,vector<T2> &v2)
{
	v2.resize(v1.size());
	v2[0]=v1[0];
	for(unsigned int i=1;i<v1.size();++i)
		v2[i]=v2[i-1]+v1[i];
}




template <class InputIterator,class InputIterator1>
void totalCumul(InputIterator first, InputIterator last,InputIterator1 first1)
{
	*first1=*first;
	while (++first != last)
	{
		*(++first1)=(*first)+(*(first1-1));
	}
}



// order brute force requires O(2N) memory
// order a vector P with id array Id
template <class T,class T1>
void orderv(vector<T> &P,vector<T1> &Id)
{
	vector<T> P1(P);
	for(unsigned int i = 0; i < P.size(); ++i)
		P[i] = P1[Id[i]];
}

// order brute force requires O(2N) memory
// order array elements PB to PE with id array Id
template <class T,class T1>
void orderf(T PB,T PE,T1 Id)
{
	int memSize=sizeof(PB[0]);
	char* P1= new char[memSize*(PE-PB)];
	memcpy(P1,&PB[0],memSize*(PE-PB));
	while(PB<PE)
	{
		memcpy(&PB[0],&P1[*(Id++)*memSize],memSize);
		PB++;
	}
	delete [] P1;
}

// order brute force requires O(N) memory
// order array elements PB to PE with id array Id
template <class T,class T1>
void order(T PB,T PE,T1 Id)
//void order(type* P,int *Id,int Nmax, int reverse)
{
	int i,Nmax=PE-PB;
	int memSize=sizeof((*PB));
	char* Psave= new char[memSize];
	int idsource, idest;
	vector<int> Id1(Nmax);
	for (i = 0; i < Nmax; i++)
		Id1[i] = Id[i];

	for (i = 0; i < Nmax; i++)
	{
		if (Id1[i] != i)
		{
			memcpy(Psave,&PB[i],memSize);
//			Psave = P[i];
			idest = i;
			while (1)
			{
				idsource = Id1[idest];
				if (idsource == i)
				{
					memcpy(&PB[idest],Psave,memSize);
//					P[idest] = Psave;
					Id1[idest] = idest;
					break;
				}
				memcpy(&PB[idest],&PB[idsource],memSize);
//				P[idest] = P[idsource];
				Id1[idest] = idest;
				idest = idsource;
			}
		}
	}
	delete [] Psave;
}


class IStats
{
public:
	IStats(){reset();}
	~IStats() {};
	void reset()
	{
		xmin=0.0;
		xmax=0.0;
		xstd=0.0;
		count=0;
	}
	void push_back(double x)
	{
		if (count> 0)
		{
			if (xmin > x)
				xmin = x;
			if (xmax < x)
				xmax = x;
			xmean += x;
			xstd += x*x;
			count++;
		}
		else
		{
			xmin=x;
			xmax=x;
			xmean=x;
			xstd=x*x;
			count=1;
		}
	}
	void print()
	{
		cout<<setw(12)<<xmin<<setw(12)<<xmax<<setw(12)<<xmean/count<<setw(12)<<sqrt(xstd*count-xmean*xmean)/count<<endl;
	}
	double xmin,xmax,xmean,xstd;
	int count;
};



#endif








