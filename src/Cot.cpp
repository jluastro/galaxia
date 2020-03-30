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

#include "Cot.h"


Cot::Cot()
{
	// TODO Auto-generated constructor stub

}

Cot::~Cot()
{
	// TODO Auto-generated destructor stub
}


void Cot::xyz_to_lbr(double *PosC,int dim)
{
	x[0]=PosC[0];	x[1]=PosC[1];	x[2]=PosC[2];
	if(dim==6)
	{
		x[3]=PosC[3];	x[4]=PosC[4];	x[5]=PosC[5];
	}
	xyz_to_lbr(&x[0],PosC,dim);
}
void Cot::xyz_to_lbr(double *PosC,double *PosS,int dim)
{
    PosS[2]=sqrt(PosC[0]*PosC[0]+PosC[1]*PosC[1]+PosC[2]*PosC[2]);
    PosS[0]=atan2(PosC[1],PosC[0]);
    if(PosS[2]>0.0)
	PosS[1]=asin(PosC[2]/PosS[2]);
    else
	PosS[1]=0.0;
    if(dim==6)
    {
    	setTM(PosS[0],PosS[1]);
    	TM.Tmult((PosC+3),(PosS+3));
    }
}
void Cot::lbr_to_xyz(double *PosS,int dim)
{
	x[0]=PosS[0];	x[1]=PosS[1];	x[2]=PosS[2];
	if(dim==6)
	{
		x[3]=PosS[3];	x[4]=PosS[4];	x[5]=PosS[5];
	}
	lbr_to_xyz(&x[0],PosS,dim);
}

// convert galactic long and latitude to Right accession and declination
// input and output units are both in degrees
void Cot::lb_to_radec(double l,double b,double &alpha,double& delta)
{
	double DTOR=3.14159/180.0;
	double RADEG=180.0/3.14159;
	double al_GP=192.85948*DTOR;
	double de_GP=27.12825*DTOR;
	double l_CP=122.932*DTOR;
	l=l*DTOR*1.0;
	b=b*DTOR*1.0;
	alpha=al_GP+atan2(cos(b)*sin(l_CP-l),cos(de_GP)*sin(b)-sin(de_GP)*cos(b)*cos(l_CP-l));
	delta=sin(alpha-al_GP)*(sin(de_GP)*sin(b)+cos(de_GP)*cos(b)*cos(l_CP-l))/(cos(b)*sin(l_CP-l));
	alpha*=RADEG;
	delta=atan(delta)*RADEG;
    if(alpha > 360.0)
    	alpha=alpha-360.0;
}

//  rad,rad  kpc 3 km/s to 3 kpc 3 km.s
// longitude l --> -Pi to +Pi ,  lat b --> -Pi/2  to Pi/2
void Cot::lbr_to_xyz(double *PosS,double *PosC,int dim)
{
    PosC[0]=PosS[2]*cos(PosS[1])*cos(PosS[0]);
    PosC[1]=PosS[2]*cos(PosS[1])*sin(PosS[0]);
    PosC[2]=PosS[2]*sin(PosS[1]);
    if(dim==6)
    {
    	setTM(PosS[0],PosS[1]);
    	TM.mult((PosS+3),(PosC+3));
    }
}



void Cot::xyz_to_lzr(double *PosC,int dim)
{
	x[0]=PosC[0];	x[1]=PosC[1];	x[2]=PosC[2];
	if(dim==6)
	{
		x[3]=PosC[3];	x[4]=PosC[4];	x[5]=PosC[5];
	}
	xyz_to_lzr(&x[0],PosC,dim);
}
void Cot::xyz_to_lzr(double *PosC,double *PosS,int dim)
{
	PosS[1]=PosC[2];
    PosS[2]=sqrt(PosC[0]*PosC[0]+PosC[1]*PosC[1]);
    PosS[0]=atan2(PosC[1],PosC[0]);
    if(dim==6)
    {
    	setTM_lzr(PosS[0],PosS[1]);
    	TM.Tmult((PosC+3),(PosS+3));
    }
}
void Cot::lzr_to_xyz(double *PosS,int dim)
{
	x[0]=PosS[0];	x[1]=PosS[1];	x[2]=PosS[2];
	if(dim==6)
	{
		x[3]=PosS[3];	x[4]=PosS[4];	x[5]=PosS[5];
	}
	lzr_to_xyz(&x[0],PosS,dim);
}
void Cot::lzr_to_xyz(double *PosS,double *PosC,int dim)
{
    PosC[0]=PosS[2]*cos(PosS[0]);
    PosC[1]=PosS[2]*sin(PosS[0]);
    PosC[2]=PosS[1];
    if(dim==6)
    {
    	setTM_lzr(PosS[0],PosS[1]);
    	TM.mult((PosS+3),(PosC+3));
    }
}





void Cot::setTM(double l1,double b1)
{
	if((l!=l1)||(b!=b1))
	{
		l=l1; b=b1;
		sin_l=sin(l);    cos_l=cos(l); sin_b=sin(b);    cos_b=cos(b);
		TM[0][0]=-sin_l;    TM[0][1]=-sin_b*cos_l;  TM[0][2]= cos_b*cos_l;
		TM[1][0]= cos_l;    TM[1][1]=-sin_b*sin_l;  TM[1][2]= cos_b*sin_l;
		TM[2][0]=   0.0;    TM[2][1]= cos_b;        TM[2][2]= sin_b;
	}
}

void Cot::setTM_lzr(double l1,double b1)
{
	if((l!=l1)||(b!=b1))
	{
		l=l1; b=b1;
		sin_l=sin(l);    cos_l=cos(l);
		TM[0][0]=-sin_l;    TM[0][1]=0.0;  TM[0][2]= cos_l;
		TM[1][0]= cos_l;    TM[1][1]=0.0;  TM[1][2]= sin_l;
		TM[2][0]=   0.0;    TM[2][1]= 1.0;        TM[2][2]= 0.0;
	}
}


const Matrix<double>& Cot::RotationMatrix(double* Pos,double theta)
{
	r=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]+Pos[2]*Pos[2]);
	x[0]=Pos[0]/r;	x[1]=Pos[0]/r;	x[2]=Pos[0]/r;
	sin_l=sin(theta);    cos_l=cos(theta);
	TR[0][0]=(1-cos_l)*x[0]*x[0]+cos_l;      TR[0][1]=(1-cos_l)*x[0]*x[1]-sin_l*x[2]; TR[0][2]=(1-cos_l)*x[0]*x[2]+sin_l*x[1];
	TR[1][0]=(1-cos_l)*x[0]*x[1]+sin_l*x[2]; TR[1][1]=(1-cos_l)*x[1]*x[1]+cos_l;      TR[1][2]=(1-cos_l)*x[1]*x[2]-sin_l*x[0];
	TR[2][0]=(1-cos_l)*x[0]*x[2]-sin_l*x[1]; TR[2][1]=(1-cos_l)*x[1]*x[2]+sin_l*x[0]; TR[2][2]=(1-cos_l)*x[2]*x[2]+cos_l;
	return TR;
}

const Matrix<double>& Cot::EulerMatrix(double phi,double theta, double psi)
{
	double cphi=cos(phi);       double sphi=sin(phi);
	double ctheta=cos(theta);   double stheta=sin(theta);
	double cpsi=cos(psi);       double spsi=sin(psi);

	TR[0][0]=cpsi*cphi-ctheta*sphi*spsi;
	TR[0][1]=cpsi*sphi+ctheta*cphi*spsi;
	TR[0][2]=spsi*stheta;

	TR[1][0]=-spsi*cphi-ctheta*sphi*cpsi;
	TR[1][1]=-spsi*sphi+ctheta*cphi*cpsi;
	TR[1][2]=cpsi*stheta;

	TR[2][0]=stheta*sphi;
	TR[2][1]=-stheta*cphi;
	TR[2][2]=ctheta;
	return TR;
}

const Matrix<double>& Cot::YPRMatrix(double phi,double theta, double psi)
{
	double cphi=cos(phi);       double sphi=sin(phi);
	double ctheta=cos(theta);   double stheta=sin(theta);
	double cpsi=cos(psi);       double spsi=sin(psi);

	TR[0][0]=ctheta*cphi;
	TR[0][1]=ctheta*sphi;
	TR[0][2]=-stheta;

	TR[1][0]=spsi*stheta*cphi-cpsi*sphi;
	TR[1][1]=spsi*stheta*sphi+cpsi*cphi;
	TR[1][2]=ctheta*spsi;

	TR[2][0]=cpsi*stheta*cphi+spsi*sphi;
	TR[2][1]=cpsi*stheta*sphi-spsi*cphi;
	TR[2][2]=ctheta*cpsi;
	return TR;
}


void Cot::rotate(double *PosC)
{
	x[0]=PosC[0];	x[1]=PosC[1];	x[2]=PosC[2];
	TR.mult(&x[0],PosC);
}

const Matrix<double>& Cot::RotationMatrix(int axis,double theta)
{
	int axis1=(axis+1)%3;
	int axis2=(axis+2)%3;
	TR.identity();
//	for(int j=0;j<3;++j)
//		for(int k=0;k<3;++k)
//			TR[j][k]=0.0;
	TR[axis1][axis1]=cos(theta);	 TR[axis1][axis2]=sin(theta);
	TR[axis2][axis1]=-sin(theta);    TR[axis2][axis2]=cos(theta);
	return TR;
}


double Cot::l=0;
double Cot::b=0;
double Cot::r=0;
double Cot::sin_b=0;
double Cot::cos_b=0;
double Cot::sin_l=0;
double Cot::cos_l=0;
std::vector<double> Cot::x=std::vector<double>(6,0.0);
Matrix<double> Cot::TM=Matrix<double>(3,3);
Matrix<double> Cot::TR=Matrix<double>(3,3);
