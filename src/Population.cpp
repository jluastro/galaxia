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

#include "Population.h"

namespace
{

double Chabrier_exp_imf(double x)
{
	return 22.8978*(exp(-pow((716.4/x),0.25)))*pow(x,-3.3);
}
double exp_imf0(double x)
{
	if(x <1)
		return pow(x,-1.6);
	else
		return pow(x,-3.0);
}




double exp_imf7(double x)
{
	return pow(x,-0.5);
}

double exp_imf8(double x)
{
	return pow(x,-0.5);
}

double exp_imf9(double x)
{
	return pow(x,-2.35);
}

//static Interp vcirc("vcirc.dat");

}

//Interp Population::vcirc=Interp("data/vcirc.dat");


void Population::setParams(int i)
{
//	double epsv[] ={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0,0.76,0.0};
//	double sigmauv[] ={16.7,19.8,27.2,30.2,36.7,43.1,43.1,67,100.0,113};
//	double sigmavv[] ={10.8,12.8,17.6,19.5,23.7,27.8,27.8,51,106,115};
//	double sigmawv[] ={6,8,10,13.2,15.8,17.4,17.5,42,85,100};
	double mass_splitv[] ={0.1e4,1.0e4,0.5e4,0.5e4,1e4,1e4,1e4,0.5e7,0.5e7,2.5e4};
	double l_splitv[] ={0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.025,0.025,0.025};
//	double den_corr[]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
//	double rho0v[]={3.09164e+06,8.75410e+06,7.36031e+06,5.00928e+06,7.93136e+06,7.03214e+06,1.04682e+07/1.25,1.34e6*1157.68,9.32e3*1407.52,0.255*13.7e9};

//	eps=epsv[i];
//	sigma_v[0]=sigmauv[i];
//	sigma_v[1]=sigmavv[i];
//	sigma_v[2]=sigmawv[i];

	mass_split=mass_splitv[i];
	l_split=l_splitv[i];

	d0=1.0;
//	rho0=rho0v[i]*den_corr[i];
//	rho_fac=den_corr[i];
	rho_fac=1.0;

	switch (i)
	{
	case 7:
		imfP= new Sampler(10000,0.07,100.0,exp_imf7,1);
//		imfP= new Sampler(10000,0.07,100.0,Chabrier_exp_imf);
		break;
	case 8:
		imfP= new Sampler(10000,0.07,100.0,exp_imf8,1);
//		imfP= new Sampler(10000,0.07,100.0,Chabrier_exp_imf,1);
		break;
	case 9:
		imfP= new Sampler(10000,0.07,100.0,exp_imf9,1);
		break;
	default:
		imfP= new Sampler(10000,0.07,100.0,exp_imf0,1);
//		imfP= new Sampler(10000,0.07,100.0,Chabrier_exp_imf);
		break;
	}

	//	cout<<i<<" "<<int_tabulated(imfP->x,imfP->px)<<" "<<imfP->meanx<<" "<<endl;
//	vcirc.print();

}

//void Population::generateVel(double* Pos)
//{
//	double r_c=radiusC(Pos);
//	double v_c=226.0;
//	v_c=vcircP->f(r_c);
//	Pos[3]=grandomn()*sigma_v[1];
//	if(ID!=8)
//		Pos[3]-=(v_c-asymmetricDrift(Pos,r_c,v_c));
//	Pos[4]=grandomn()*sigma_v[2];
//	Pos[5]=grandomn()*sigma_v[0];
//	Cot::xyz_to_lzr(Pos);
//	Cot::lzr_to_xyz(Pos,6);
////	Pos[4]-=226.0;
//}

void Population::setVel(StarParticle &Star)
{
	double age_temp=pow(10.0,Star.age());
	double r_c=radiusC(&Star.pos(0));
	double v_c=226.0;
	v_c=vcircP->f(r_c);
	double temp=sigmar_func(age_temp,r_c);

	Star.pos(3)=grandomn()*temp*sigma_v[1]/sigma_v[0]-(v_c-asymmetricDrift(&Star.pos(0),v_c,age_temp,temp)*1.0);
	Star.pos(4)=grandomn()*temp*sigma_v[2]/sigma_v[0];
	Star.pos(5)=grandomn()*temp;

//	if(ID<=7)
//		Star.Pos[3]=-(v_c+vphi_dist(age_temp,r_c));

	if(ID==8)
	{
		Cot::xyz_to_lbr(&Star.pos(0));
		Cot::lbr_to_xyz(&Star.pos(0),6);
	}
	else
	{
		Cot::xyz_to_lzr(&Star.pos(0));
		Cot::lzr_to_xyz(&Star.pos(0),6);
	}

//	cout<<r_c<<" "<<v_c<<" "<<Star.Pos[4]-v_c<<endl;
//	Pos[4]-=226.0;
}



Population::Population()
{
	Rsun=8.0;
//	PosC[0]=-Rsun; PosC[1]=0.0 ; PosC[2]=0.015;
	PosC[0]=-Rsun; PosC[1]=0.0 ; PosC[2]=0.015;
	// TODO Auto-generated constructor stub
	vconst=226.84;
}

Population::~Population()
{
	if(imfP!=NULL)
		delete imfP;
}

ExpDisk::~ExpDisk()
{
}

ThinDisk::~ThinDisk()
{
}


ThinDisk0::~ThinDisk0()
{
}

ThinDisk1::~ThinDisk1()
{
}

ThickDisk::~ThickDisk()
{
}

Spheroid::~Spheroid()
{
}

Bulge::~Bulge()
{
}
