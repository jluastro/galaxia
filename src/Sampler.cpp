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



//#include<cmath>
//#include<iostream>
//#include<fstream>
#include "Sampler.h"
#include "Functions.h"
//#include "psplot.h"
//using namespace std;


//Sampler::Sampler(int nsize,double xmin1,double xmax1,double (* func) (double x))
//{
////	 srand48(4);
//	 function=func;
////	x.resize(nsize);
////	 linspace(x,xmin1,xmax1,nsize);
//	loglinspace(x,xmin1,xmax1,nsize);
//	px.resize(x.size());
//	for(int i=0;i<int(x.size());++i)
//		px[i]=function(x[i]);
//	normalize();
//}

Sampler::Sampler(int nsize,double xmin1,double xmax1,double (* func) (double x),int optionlinlog)
{
//	 srand48(4);
	 function=func;
//	x.resize(nsize);
	 if(optionlinlog==0)
		 linspace(x,xmin1,xmax1,nsize);
	 else
	 {
		 certify(xmax1>xmin1,"Smapler test failed xmin1>xmax1");
		 if(xmin1<=0)
		 {
			 cout<<"Sampler:log spacing not allowed for xmin<=0"<<endl;
			 exit(1);
		 }
//		 loglinspace(x,xmin1,xmax1,nsize);
		 logspace(x,xmin1,xmax1,nsize);
	 }
//	 printv(x.begin(),x.end(),"Sampler");
	px.resize(x.size());
	for(int i=0;i<int(x.size());++i)
		px[i]=function(x[i]);
	normalize();
}

void Sampler::print( )
{
	cout<<"Printing Sampler Data "<<endl;
	for(unsigned int i=0;i<x.size();++i)
		cout<<i<<" x="<<x[i]<<" px="<<px[i]<<" cpd="<<cpd[i]<<endl;
}

void Sampler::plot( )
{
//	PSplot plot1;
//	plot1.plot(x,px,"x","p(x)",1,1,0);
//	plot1.oplot(x,cpd,2,0.5);
//	plot1.close();
}

Sampler::Sampler(const string fname)
{
	   float temp;
	    ifstream fd;
	    int i,no;

	    fd.open(fname.c_str());
	    if (fd.is_open())
	    {
	    	cout<<"Reading luminosity function file: "<<fname<<endl;
		fd>>no;
		for(i=0;i<no;++i)
		{
		    fd>>temp;
		    x.push_back(temp);
		}
		for(i=0;i<no;++i)
		{
		    fd>>temp;
		    px.push_back(temp);
		}
		fd>>temp;

		fd.close();
	    }
	    else
		cout<<"Error opening file"<<endl;
		   cout<<"Int P(l)dl "<<int_tabulated(x,px)<<endl;
		normalize();
}

Sampler::Sampler(vector<double> &x1,vector<double> &px1)
{
	x=x1;
	px=px1;
	normalize();
	calculateCpd();
}

Sampler::Sampler(vector<double> &x1)
{
	x=x1;
	sort(x.begin(),x.end());
	double temp;
	cpd.resize(x.size());
	for(unsigned int i=0;i<x.size();++i)
    {
	    cpd[i]=i*1.0/(x.size()-1);
    }
	vector<double>::iterator it=unique2(x.begin(),x.end(),cpd.begin());
	cpd.resize(it-x.begin());
	x.resize(it-x.begin());
	cpd.back()=1.0;
	px.resize(x.size());
	for(unsigned int i=1;i<(x.size()-1);++i)
    {
	    temp=(cpd[i]-cpd[i-1])/(x[i]-x[i-1])+(cpd[i+1]-cpd[i])/(x[i+1]-x[i]);
	    temp*=0.5;
	    px[i]=temp;
    }
    px[0]=0.0;
    px[x.size()-1]=0.0;
	x_min=*min_element(x.begin(),x.end());
	x_max=*max_element(x.begin(),x.end());
}


Sampler::~Sampler()
{
}



//void Sampler::linspace(vector<double>&x,double xmin,double xmax)
//{
//	double temp=(xmax-xmin)/(x.size()-1);
//	x.front()=xmin;
//	for(unsigned int i=1;i<x.size();++i)
//		x[i]=xmin+i*temp;
//
//}
//
//
//
//
//void Sampler::logspace(vector<double> &x,double xmin,double xmax)
//{
//	double temp=(log10(xmax)-log(xmin))/(x.size()-1);
//	x.front()=xmin;
//	for(unsigned int i=1;i<x.size();++i)
//		x[i]=pow(10.0,log10(xmin)+i*temp);
//}

//double Sampler::int_tabulated(vector<double> &x,vector<double> &y,double xmin,double xmax)
//{
//    double temp=0;
//    for(unsigned int i=1;i<x.size();++i)
//    {
//    	if(((x[i])<=xmax)&&(x[i-1]>=xmin))
//    		temp+=(x[i]-x[i-1])*(y[i]+y[i-1])*0.5;
//    }
//    return temp;
//}

//double Sampler::int_tabulated(vector<double> &x,vector<double> &y)
//{
//	double temp=0;
//	unsigned int i;
//    for(i=1;i<x.size();++i)
//    {
//	    temp+=(x[i]-x[i-1])*(y[i]+y[i-1])*0.5;
//    }
////    for(i=4;i<x.size();i+=4)
////    {
////	    temp+=(x[i]-x[i-4])*(7*y[i-4]+32*y[i-3]+12*y[i-2]+32*y[i-1]+7*y[i])*4.0/90.0;
////    }
////    for(unsigned int j=i-4;j<x.size();++j)
////    {
////	    temp+=(x[j]-x[j-1])*(y[j]+y[j-1])*0.5;
////    }
//
//    return temp;
//}

void Sampler::normalize( )
{
	double norm=int_tabulated(x,px)	;
//	cout<<"Normalization factor "<<norm<<" "<<x.size()<<endl;
//	norm=1.0;
	for(unsigned int i=0;i<x.size();++i)
	px[i]/=norm;
	meanx=0.0;
	vector<double> temp(x.size(),0.0);
	for(unsigned int i=0;i<x.size();++i)
		temp[i]=px[i]*x[i];
	meanx=int_tabulated(x,temp);
	x_min=*min_element(x.begin(),x.end());
	x_max=*max_element(x.begin(),x.end());
	calculateCpd();
	//	cout<<"Mean x "<<meanx<<endl;
}




void Sampler::setSeed(int64_t seed)
{
	srand48(seed)	;
}

void Sampler::setRange(double xmin1,double xmax1)
{
	unsigned int	up=0;
	if(xmax1 < xmin1)
	{
		double temp=xmin1;
		xmin1=xmax1;
		xmax1=temp;
	}
	if(xmax1>=x_max)
	{
		cpd_max=cpd.back();
	}
	else	if(xmax1<=x_min)
	{
		cpd_max=cpd.front();
	}
	else
	{
		up=(upper_bound(x.begin(),x.end(),xmax1)-x.begin())-1;
		cpd_max=cpd[up]+(cpd[up+1]-cpd[up])*(xmax1-x[up])/(x[up+1]-x[up]);
	}


	if(xmin1>=x_max)
	{
		cpd_min=cpd.back();
	}
	else	if(xmin1<=x_min)
	{
		cpd_min=cpd.front();
	}
	else
	{
		up=(upper_bound(x.begin(),x.end(),xmin1)-x.begin())-1;
		cpd_min=cpd[up]+(cpd[up+1]-cpd[up])*(xmin1-x[up])/(x[up+1]-x[up]);
	}

// change done Jul 28 ,2009
	if(xmax1<=x_min)
		cpd_min=cpd[1];
	if(xmin1>=x_max)
		cpd_max=cpd[cpd.size()-2];


	if(cpd_min==cpd_max)
	{
//		cpd_max=cpd[up+1];
		cout<<"Warning:Range may be too small "<<xmin1<<" "<<xmax1<<" "<<x.front()<<" "<<x.back()<<endl;
		exit(1);
	}

//	range_fac=0.0;
//	int lo=int(find_if(x.begin(),x.end(),bind2nd(greater_equal<double>(),xmin1))-x.begin());
//	up=int(find_if(x.begin(),x.end(),bind2nd(greater<double>(),xmax1))-x.begin());
//	cout<<xmin1<<" "<<lo<<" "<<x[lo]<<" "<<xmax1<<" "<<up<<" "<<x[up]<<endl;
//	if( lo==0) lo=1;
//	for(int i=lo;i<up;++i)
//		range_fac+=(px[i]+px[i-1])*(x[i]-x[i-1])*0.5;
//	cout<<"Range fac "<<range_fac<<endl;

}

void Sampler::getFacv(vector<double> &x1_a,vector<double> &x2_a,vector<double> &fac_a)
{
//	double range_fac=0.0;
	fac_a.resize(x1_a.size());
//	fac_a[0]=0.0;
//	cout<<x.front()<<x.back()<<" "<<cpd.front()<<" "<<cpd.back()<<endl;
	for(unsigned int j=0;j<(x1_a.size());++j)
		fac_a[j]=fabs(interpolate(x,cpd,x2_a[j])-interpolate(x,cpd,x1_a[j]));
}



//void Sampler::getFacv(vector<double> &x1_a,vector<double> &x2_a,vector<double> &fac_a)
//{
//	double range_fac=0.0;
//	fac_a.resize(x1_a.size());
//	fac_a[0]=0.0;
//	for(unsigned int j=1;j<(x1_a.size());++j)
//	{
//		range_fac=0.0;
//		int lo=int(find_if(x.begin(),x.end(),bind2nd(greater_equal<double>(),min(x1_a[j],x1_a[j-1])))-x.begin());
//		int up=int(find_if(x.begin(),x.end(),bind2nd(greater<double>(),max(x2_a[j],x2_a[j-1])))-x.begin());
////		cout<<j<<" "<<lo<<" "<<up<<" "<<x.size()<<endl;
//		if(up>(x.size()-1))
//			up=x.size()-1;
//		if(lo>(x.size()-1))
//			lo=x.size()-1;
//		if(lo<=0)
//			lo=1;
//
//		for(int i=lo;i<=up;++i)
//			range_fac+=(px[i]+px[i-1])*(x[i]-x[i-1])*0.5;
//		fac_a[j]=range_fac;
////		cout<<j<<" here  "<<range_fac<<endl;
//	}
////	cout<<"Range fac "<<range_fac<<endl;
//
//}

double Sampler::getFac(double xmin1,double xmax1)
{
	return fabs(interpolate(x,cpd,xmin1)-interpolate(x,cpd,xmax1));
}


double Sampler::rand( )
{
	double u=cpd_min+grandomu()*(cpd_max-cpd_min);
	int up;
	if(cpd.front()==0.0)
		up=int(upper_bound(cpd.begin(),cpd.end(),u)-cpd.begin())-1;
	else
		up=int(upper_bound(cpd.begin(),cpd.end(),u,greater<double>())-cpd.begin())-1;
	//	cout<<cpd.front()<<" a "<<cpd.back()<<" "<<up<<" "<<cpd[up+1]<<" "<<cpd[up]<<endl;
	return x[up]+(x[up+1]-x[up])*(u-cpd[up])/(cpd[up+1]-cpd[up]);
}

vector<double> Sampler::randv(int nsize1)
{
	vector<double> y(nsize1);
	for(int i=0;i<nsize1;++i)
		y[i]=rand();
	return y;
}

void Sampler::calculateCpd( )
{
	cpd.resize(px.size());
	if(px.front()<px.back())
	{
		cpd.front()=0.0;
		double temp=0.0;
		for(unsigned int i=1;i<x.size();++i)
		{
			temp+=(x[i]-x[i-1])*(px[i]+px[i-1])*0.5;
			cpd[i]=temp;
//			cout<<i<<" "<<cpd[i]<<" "<<(cpd[i]-cpd[i-1])<<endl;
		}
	}
	else
	{
		cpd.back()=0.0;
		double temp=0.0;
		for(int i=x.size()-2;i>=0;--i)
		{
			temp+=(x[i+1]-x[i])*(px[i+1]+px[i])*0.5;
			cpd[i]=temp;
//			cout<<i<<" "<<cpd[i]<<" "<<(cpd[i]-cpd[i+1])<<endl;
		}
	}

	cpd_min=cpd.front();
	cpd_max=cpd.back();
}











