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

#ifndef GEOMETRY_H_
#define GEOMETRY_H_
#include "StarParticle.h"


double periodic(double dx,double l);

class Geometry
{
public:
	Geometry();
	virtual int checkSphere(const double* Pos,const double &h)=0;
	virtual int checkP(const double* Pos)=0;
	virtual ~Geometry();
protected:
};

class AllSky:public Geometry
{
public:
	AllSky();
	~AllSky();
	int checkSphere(const double* Pos,const double &h){return 1;}
	int checkP(const double* Pos) {return 1;}
};



class Cone:public Geometry
{
public:
	Cone(double b_Deg,double db_Deg)
	{
		if((b_Deg)<0.0)
		{
			cout<<"Range less than 0.0"<<endl;
			exit(1);
		}
		if((b_Deg+db_Deg)>180.0)
		{cout<<"Range greater than 180.0"<<endl; exit(1);}
		th0=max((b_Deg)*PI/180.0,0.0);
		th1=min((b_Deg+db_Deg)*PI/180.0,PI);
	}
	int checkSphere(const double* Pos,const double &h);
	int checkP(const double* Pos);
	~Cone();
private:
//	double PosC[3];
	double th0,th1;
};



class Wedge:public Geometry
{
public:
	Wedge(double b_Deg,double db_Deg,int axis1)
	{
		axis=axis1;
		if((b_Deg)<0.0)
		{
			cout<<"Range less than 0.0"<<endl;
			exit(1);
		}
		if((b_Deg+db_Deg)>180.0)
		{cout<<"Range greater than 180.0"<<endl; exit(1);}
		th0=max((b_Deg)*PI/180.0,0.0);
		th1=min((b_Deg+db_Deg)*PI/180.0,PI);
	}
	int checkSphere(const double* Pos,const double &h);
	int checkP(const double* Pos);
	~Wedge();
private:
//	double PosC[3];
	double th0,th1;
	int axis;
};



class Strip:public Geometry
{
public:
	Strip(double l_Deg,double dl_Deg)
	{
		if((l_Deg)<0.0)
		{cout<<"Range less than 0.0"<<endl; exit(1);}
		if((l_Deg)>360.0)
		{cout<<"Range greater than 360.0"<<endl; exit(1);}
		ph0=(l_Deg)*PI/180.0;
		ph1=(l_Deg+dl_Deg)*PI/180.0;
	}
	double periodic(double dx,double l);
	int checkSphere(const double* Pos,const double &h);
	int checkP(const double* Pos);
	~Strip();
private:
	double ph0,ph1;
};


class Patch:public Geometry
{
public:
	Patch(double omegaSqDeg)
	{
		if(omegaSqDeg<0)
		{
			th0=0.0;
			th1=PI;
		}
		else
		{
			th0=0.0;
			if((omegaSqDeg*(PI*PI)/(2*PI*180.0*180.0))>2)
				th1=PI;
			else
				th1=acos(1-omegaSqDeg*(PI*PI)/(2*PI*180.0*180.0));
		}
	}
	int checkSphere(const double* Pos,const double &h);
	int checkP(const double* Pos) {return 1;}
	~Patch();
private:
//	double PosC[3];
	double th0,th1;
};


class GPatch:public Geometry
{
public:
	GPatch(double l,double b,double d_th1)
	{
		PosU[0]=cos(b*PI/180.0)*cos(l*PI/180.0);
		PosU[1]=cos(b*PI/180.0)*sin(l*PI/180.0);
		PosU[2]=sin(b*PI/180.0);
		d_th=d_th1*PI/180.0;
	}
	int checkSphere(const double* Pos,const double &h)
	{
		tempr = sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1]+ Pos[2] * Pos[2]);
		if (h >= tempr)
			return 1;
		temp=(Pos[0] * PosU[0] + Pos[1] * PosU[1] + Pos[2] * PosU[2])/tempr;
		if((acos(temp)-asin(h/tempr))<d_th)
			return 1;
		else
			return 0;

	}
	int checkP(const double* Pos)
	{
		tempr = sqrt(Pos[0] * Pos[0] + Pos[1] * Pos[1] + Pos[2] * Pos[2]);
		temp=(Pos[0] * PosU[0] + Pos[1] * PosU[1] + Pos[2] * PosU[2])/tempr;
		if(acos(temp)<d_th)
			return 1;
		else
			return 0;
	}
	~GPatch();
private:
	double PosU[3];
	double d_th,tempr,temp;
};



#endif /* GEOMETRY_H_ */
