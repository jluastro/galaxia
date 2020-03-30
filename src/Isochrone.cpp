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

#include <functional>
#include <algorithm>
#include "Isochrone.h"
#include "Functions.h"

Isochrone::Isochrone()
{
}

Isochrone::~Isochrone()
{
}

void Isochrone:: print()
{
	cout<<"Age "<<age<<" FeH "<<FeH<<" Alpha "<<alpha<<endl;
	for(size_t i=0;i<m.size();++i)
		cout<<setw(16)<<m[i]<<setw(16)<<Mag0[i]<<setw(16)<<Mag1[i]<<endl;
}

void Isochrone:: setMon()
{
	Vmon=Mag0;
	// set monotonic decreasing abs mag array Vmon
	for(unsigned int i=1;i<m.size();++i)
	{
		if(Vmon[i]>Vmon[i-1])
			Vmon[i]=Vmon[i-1];
	}

	maxV=Mag0.front();
	m_min=m.front();
	m_max=m.back();
}

void Isochrone:: setTip()
{
	mtip=m.front();
	for(unsigned int i=1;i<m.size();++i)
	{
		if((Lum[i]<Lum[i-1])&&(Teff[i]>Teff[i-1])&&(m[i]>0.15))
		{
			mtip=m[i];
			break;
		}
	}
}


void Isochrone:: setmaxV(double maxV1)
{
	maxV=maxV1;
	if(maxV>=Vmon.front())
		m_min=m.front();
	else if(maxV<=Vmon.back())
		m_min=m.back();
	else
	{
		int i = int(upper_bound(Vmon.begin(), Vmon.end(), maxV, greater<double> ())- Vmon.begin() - 1);
		m_min=m[i];
	}
}

void Isochrone:: interpolate(StarParticle &Star)
{
	float smass_t=Star.smass();
//	if(Star.smass()<=m.front())
	if(smass_t<=m_min)
	{
		Star.mag(0)=1e6;
		Star.mag(1)=1e6;
		Star.mag(2)=1e6;
		Star.lum()=1e6;
		Star.teff()=1e6;
		Star.grav()=1e6;
	}
	else if(smass_t>=m.back())
	{
			Star.mag(0)=1e6;
			Star.mag(1)=1e6;
			Star.mag(2)=1e6;
			Star.lum()=1e6;
			Star.teff()=1e6;
			Star.grav()=1e6;
	}
	else
	{
		int i=locate(m,double(smass_t));
		double temp=(smass_t-m[i])/(m[i+1]-m[i]);
		Star.mag(0)=Mag0[i]+(Mag0[i+1]-Mag0[i])*temp;
		Star.mag(1)=Mag1[i]+(Mag1[i+1]-Mag1[i])*temp;
		Star.mag(2)=Mag2[i]+(Mag2[i+1]-Mag2[i])*temp;
		Star.lum()=Lum[i]+(Lum[i+1]-Lum[i])*temp;
		Star.teff()=Teff[i]+(Teff[i+1]-Teff[i])*temp;
		Star.grav()=Grav[i]+(Grav[i+1]-Grav[i])*temp;
	}
//	cout<<Star.mag(0)<<" "<<Star.mag(1)<<endl;
}

void Isochrone::interpolateTGM(float smass,vector<double> &x)
{
	if(x.size()>(Mags.size()+5))
	{
		cout<<"size mismatch of x and Mags in interpolateTGM"<<endl;
	}

	x[4]=mtip;
	if(smass>=m.back())
	{
		x[0]=Lum.back();
		x[1]=Teff.back();
		x[2]=Grav.back();
		x[3]=Mact.back();
		if(x.size()>5)
			for(size_t k=5;k<x.size();++k)
				x[k]=Mags[k-5].back();
		//				x[k]=1e4;
	}
	else if(smass<=m.front())
	{
		x[0]=Lum.front();
		x[1]=Teff.front();
		x[2]=Grav.front();
		x[3]=Mact.front();
		if(x.size()>5)
			for(size_t k=5;k<x.size();++k)
				x[k]=Mags[k-5].front();
//				x[k]=1e4;
	}
	else
	{
		int i=locate(m,double(smass));
		double temp=(smass-m[i])/(m[i+1]-m[i]);
		x[0]=Lum[i]+(Lum[i+1]-Lum[i])*temp;
		x[1]=Teff[i]+(Teff[i+1]-Teff[i])*temp;
		x[2]=Grav[i]+(Grav[i+1]-Grav[i])*temp;
		x[3]=Mact[i]+(Mact[i+1]-Mact[i])*temp;

		if(x.size()>5)
		for(size_t k=5;k<x.size();++k)
		{
			x[k]=Mags[k-5][i]+(Mags[k-5][i+1]-Mags[k-5][i])*temp;
		}
	}

}


//void Isochrone:: setTip1()
//{
////	tip=m.size()-1;
////	print();
//	cout<<"Age "<<age<<" FeH "<<FeH<<" Alpha "<<alpha<<" "<<tip<<" "<<m.size()<<endl;
//	unsigned int i;
//	for(i=0;i<(m.size()-1);++i)
//	{
////		cout<<i<<" "<<m[i]<<" "<<m[i+1]<<" "<<(m[i]==m[i+1])<<endl;
//		if( (m[i]==m[i+1]) && (((B[tip+1]-V[tip+1])-(B[tip]-V[tip]))<0.0) )
//		{
//			break;
//		}
//	}
//
//	if(i==(m.size()-1))
//			cout<<"tip not found "<<tip<<" "<<m.size()-1<<endl;
//	else
//		cout<<"found "<<tip<<" "<<m.size()-1<<" "<<V[i+1]-V[i]<<" "<<(B[i+1]-V[i+1])-(B[i]-V[i])<<endl;
//
//}
