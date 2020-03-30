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

#include "Geometry.h"

Geometry::Geometry()
{
}


Geometry::~Geometry()
{
	// TODO Auto-generated destructor stub
}

//----------------------------------------

AllSky::AllSky()
{

}
AllSky::~AllSky()
{
	// TODO Auto-generated destructor stub
}

//---------------------------------------------------------
int Cone::checkSphere(const double* Pos,const double &h)
{
	double r = sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1]
			+ Pos[2] * Pos[2]);
//	cout<<h<<" "<<Part->h(1)<<" "<<r<<" "<<endl;

	if (h >= r)
		return 1;

	double th = acos(Pos[2] / r);
	double dth = asin(h / r);
//	cout<<h<<" "<<th<<" "<<dth<<" "<<th0<<" "<<th1<<endl;
	if ((th + dth) < th0)
		return 0;
	if ((th - dth) > th1)
		return 0;
	return 1;
}


int Cone::checkP(const double* Pos)
{
	double r = sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1]
			+ Pos[2] * Pos[2]);
	double th = acos(Pos[2] / r);
//	cout<<h<<" "<<th*180.0/3.14<<" "<<" "<<th0*180.0/3.14<<" "<<th1*180.0/3.14<<endl;
	if (th < th0)
		return 0;
	if (th > th1)
		return 0;
//	cout<<h<<" "<<th*180.0/3.14<<" "<<" "<<th0*180.0/3.14<<" "<<th1*180.0/3.14<<endl;
	return 1;
}
Cone::~Cone()
{
	// TODO Auto-generated destructor stub
}


//---------------------------------------------------------
int Wedge::checkSphere(const double* Pos,const double &h)
{
	double r = sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1]
			+ Pos[2] * Pos[2]);
	if (h >= r)
		return 1;

	// theta+d_theta>60.0
	if((acos(Pos[2] / r)-asin(h / r))>PI/3)
		return 0;

	double th = acos(Pos[axis] / r);
	double dth = asin(h / r);
//	cout<<h<<" "<<th<<" "<<dth<<" "<<th0<<" "<<th1<<endl;
	if ((th + dth) < th0)
		return 0;
	if ((th - dth) > th1)
		return 0;
	return 1;
}

int Wedge::checkP(const double* Pos)
{

	double r = sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1]
			+ Pos[2] * Pos[2]);
	// theta>60.0
		if((acos(Pos[2] / r))>PI/3)
			return 0;


	double th = acos(Pos[axis] / r);
//	cout<<h<<" "<<th*180.0/3.14<<" "<<" "<<th0*180.0/3.14<<" "<<th1*180.0/3.14<<endl;
	if (th < th0)
		return 0;
	if (th > th1)
		return 0;
//	cout<<h<<" "<<th*180.0/3.14<<" "<<" "<<th0*180.0/3.14<<" "<<th1*180.0/3.14<<endl;
	return 1;
}
Wedge::~Wedge()
{
	// TODO Auto-generated destructor stub
}

//---------------------------------------------------------
double Strip::periodic(double dx,double l)
{
	return fabs(dx)<(0.5*l) ? dx: ( (dx<0) ? (l+dx): (dx-l) );
//	if(fabs(dx)>(0.5*l))
//	return ( (dx<0) ? (l+dx): (dx-l) );
//	if(dx<(-l*0.5))
//		return l+dx;
//	if(dx>(l*0.5))
//		return dx-l;
//	return dx;
}
int Strip::checkSphere(const double* Pos,const double &h)
{
	double r = sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1] + Pos[2] * Pos[2]);
	if (h >= r)
		return 1;
	double th = acos(Pos[2] / r);
	double dth = asin(h / r);
//	cout<<h<<" "<<th<<" "<<dth<<" "<<th0<<" "<<th1<<endl;
	if ((th - dth) > PI)
		return 0;


	r=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]);
	if(h>=r)
		return 1;
	double ph=acos(Pos[0]/r);
	if(Pos[1]<0)
		ph=2*PI-ph;
	double dph=asin(h/r);
	if(periodic(ph+dph-ph0,2*PI)<0)
		return 0;
	if(periodic(ph-dph-ph1,2*PI)>0)
		return 0;
	return 1;
}

int Strip::checkP(const double* Pos)
{
	double r=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]);
	double ph=acos(Pos[0]/r);
	if(Pos[1]<0)
		ph=2*PI-ph;
	if(periodic(ph-ph0,2*PI)<0)
		return 0;
	if(periodic(ph-ph1,2*PI)>0)
		return 0;
	return 1;
}

Strip::~Strip()
{
	// TODO Auto-generated destructor stub
}


GPatch::~GPatch()
{
	// TODO Auto-generated destructor stub
}
