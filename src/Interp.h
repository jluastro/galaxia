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

#ifndef INTERP_H_
#define INTERP_H_


#include<algorithm>
#include "sutils.h"

template<class T>
void interpolate(vector<T>&x,vector<T> &y,vector<T>&x1,vector<T> &y1)
{
	int p;
	y1.resize(x1.size());
	for(unsigned int i=0;i<x1.size();++i)
	{
		if(x1[i]>=x.back())
			y1[i]=y.back();
		else 		if(x1[i]<x.front())
			y1[i]=y.front();
		else
		{
			p=int(upper_bound(x.begin(),x.end(),x1[i])-x.begin()-1);
			y1[i]=y[p]+(y[p+1]-y[p])*(x1[i] -x[p])/(x[p+1]-x[p]);
		}
//		cout<<i<<" "<<x1[i]<<" "<<y1[i]<<" "<<x.front()<<" "<<x.back()<<" "<<y.front()<<" "<<y.back()<<endl;
	}
}

template<class T>
int inline locate(vector<T>&x,T x1)
{
	if(x1<=x.front())
		return 0;
	if(x1>=x.back())
		return int(x.size()-2);
	return int(upper_bound(x.begin(),x.end(),x1)-x.begin()-1);
}

template<class T>
int inline locate_nearest(vector<T>&x,T x1)
{
	if(x1<=x.front())
		return 0;
	if(x1>=x.back())
		return int(x.size()-1);
	int i=int(upper_bound(x.begin(),x.end(),x1)-x.begin()-1);

	if((i<0)||(i>=(int(x.size())-1)))
	{
		cout<<i<<" "<<x.size()<<" "<<x1<<" "<<x.front()<<" "<<x.back()<<endl;
	}
	assert((i>=0)&&(i<(int(x.size())-1)));

	if (((x1-x[i])/(x[i+1]-x[i]))<0.5)
		return i;
	else
		return (i+1);
}



template<class T>
void inline interpolate(vector<T>&x,vector<T> &y, T &x1,T &y1)
{
	int p;
	if(x1>=x.back())
		y1=y.back();
	else
		if(x1<x.front())
			y1=y.front();
		else
		{
			p=int(upper_bound(x.begin(),x.end(),x1)-x.begin()-1);
			assert((p>=0)&&(p<int(x.size()-1)));
			y1=y[p]+(y[p+1]-y[p])*(x1 -x[p])/(x[p+1]-x[p]);
		}
}

//template<class T>
//T interpolate(vector<T>&x,vector<T> &y, T x1)
//{
//
//	int p;
//	if(x1>=x.back())
//	{
//		return y.back();
//	}
//	if(x1<x.front())
//	{
//		return y.front();
//	}
//	p=int(upper_bound(x.begin(),x.end(),x1)-x.begin()-1);
//	return y[p]+(y[p+1]-y[p])*(x1 -x[p])/(x[p+1]-x[p]);
//}

template<class T>
T interpolate(vector<T>&x,vector<T> &y, T x1)
{

	int stat_p=0;

	if(x1>=x.back())
	{
		return y.back();
	}
	if(x1<x.front())
	{
		return y.front();
	}

	stat_p=int(upper_bound(x.begin(),x.end(),x1)-x.begin()-1);


//	if((stat_p<0)||(stat_p>=(x.size()))||(x1<x[stat_p])||(x1>=x[stat_p+1]))
//	{
//		stat_p=int(upper_bound(x.begin(),x.end(),x1)-x.begin()-1);
//	}
//	if((stat_p<0)||(stat_p>=(x.size()-1)))
//	{
//		cout<<stat_p<<" "<<x.size()<<" "<<x1<<" "<<x.front()<<" "<<x.back()<<endl;
//	}
	assert((stat_p>=0)&&(stat_p<int(x.size()-1)));

//	{
//		cout<<stat_x<<" "<<stat_p<<" "<<x1<<" "<<x[stat_p]<<" "<<x[stat_p+1]<<" "<<y[stat_p]<<" "<<y[stat_p+1]<<endl;
//	}
		return y[stat_p]+(y[stat_p+1]-y[stat_p])*(x1 -x[stat_p])/(x[stat_p+1]-x[stat_p]);
}


template<class T>
void printv1(T a,T b)
{
	std::cout<<"Printing Vector"<<endl;
	while(a!=b)
//		std::cout<<*a++<<" "<<endl;
		std::cout<<scientific<<*a++<<endl;
//		std::cout<<scientific<<(*a)<<setw(12)<<endl;
	std::cout<<endl;
}


class Interp
{
public:
	Interp(vector<double> &x1,vector<double> &y1)
	{
		x=x1;
		y=y1;
	}
	void print()
	{
		printv1(x.begin(),x.end());
		printv1(y.begin(),y.end());
	}
	Interp(const string& fname)
	{
		setFromFile(fname);
//		cout<<"Init interp"<<endl;
	}
	Interp & operator=(const Interp &a)
	{
		x=a.x;
		y=a.y;
		return *this;
	}
	Interp();
	void setFromFile(const string& fname);
	inline double f(double x1)
	{
		return interpolate(x,y,x1);
	}
	virtual ~Interp();
	vector<double> x,y;
};

#endif /* INTERP_H_ */
