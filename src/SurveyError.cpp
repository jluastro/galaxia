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

//#include <cmath>
//#include "Functions.h"
#include "SurveyError.h"
#include "Cot.h"

void SurveyError::initialize()
{
//	TM.resize(3,3);
	appMag.resize(6);
	vr.resize(6);
	appMag[0]=5.0;  vr[0]=0.0;
	appMag[1]=10.0; vr[1]=0.1;
	appMag[2]=16.0; vr[2]=1.0;
	appMag[3]=17.0; vr[3]=2.0;
	appMag[4]=18.0; vr[4]=6.0;
	appMag[5]=20.0; vr[5]=10.0;
	sigma_w=1000.0;
	sigma_mu=200.0;
	sigma_r=0.1;
	sigma_vr=10.0;
	sigma_vlb=10.0;
	sigma_fe=0.1;
	sigma_al=0.1;
	// TODO Auto-generated constructor stub
}


void SurveyError::add(StarParticle &Star)
{
//	Star.feh()=log10(Star.feh()/0.019);
	if (errorOption > 0)
	{
		double PosG[6];
		Cot::xyz_to_lbr(&Star.pos(0), PosG);

		if (errorOption == 1)
		{
			double z=-1.0,z2=1.0;
			sigma_vlb = PosG[2] * 1e-3 * 4.7407 * sigma_mu;
			if (PosG[2] != 0.0)
			{
				while(z < 0.0)
					z= PosG[2]*(1 + sigma_r * gauss.dev());
			    z2=z/PosG[2];
				PosG[2] =z;
			}
			PosG[3] += sigma_vlb * gauss.dev(); PosG[3] *=z2;
			PosG[4] += sigma_vlb * gauss.dev(); PosG[4] *=z2;
			PosG[5] += sigma_vr * gauss.dev();
		}
		else
		if (errorOption == 2)
		{
			double V_I = Star.mag(1) - Star.mag(2);
			double G = Star.mag(0) + 0.51 - 0.5 * sqrt(0.6 + (V_I - 0.6) * (V_I - 0.6))
					- 0.065 * (V_I - 0.6) * (V_I - 0.6);
			double z = pow(10.0, 0.4 * (G - 15));
			double z2 = z * z;
			sigma_w = sqrt(7.0 + 105.0 * z + 1.3 * z2 + z2 * z2 * z2 * 6.0e-10)
					* (0.96 + 0.04 * V_I);
			sigma_vr = interpolate(appMag, vr, Star.mag(0));
			sigma_mu=0.75*sigma_w;
			sigma_vlb = PosG[2] * 1e-3 * 4.7407 * sigma_mu;
			double sigma_wlb = 0.87 * sigma_w * 1e-6 * (3.14159 / (180.0 * 3600.0));
			PosG[0] += sigma_wlb * gauss.dev();
			PosG[1] += sigma_wlb * gauss.dev();
			z2=1.0;
			if (PosG[2] != 0.0)
			{
			    if ((sigma_w*PosG[2]/1000.0)>sigma_r)
			    	z= PosG[2]*(1 + sigma_r * gauss.dev());
			    else
			    {
				z = 1000.0 / (1000.0 / PosG[2] + sigma_w* gauss.dev());
			    }
			    while(z < 0.0)
				 z= PosG[2]*(1 + sigma_r * gauss.dev());
			    z2=z/PosG[2];
			    PosG[2]=z;
			}

			PosG[3] += sigma_vlb * gauss.dev(); PosG[3] *=z2;
			PosG[4] += sigma_vlb * gauss.dev(); PosG[4] *=z2;
			PosG[5] += sigma_vr * gauss.dev();
		}
		else
		if (errorOption == 3)
		{
			double z=-1.0,z2=1.0;
			sigma_vlb = PosG[2] * 1e-3 * 4.7407 * sigma_mu;
			if (PosG[2] != 0.0)
			{
				z= pow(10.0,(Star.mag(0)-4.7-10.0)/5.0);
				if(z < 0.0)
					z=0.1;
				z2=z/PosG[2];
				PosG[2] =z;
			}
			PosG[3] += sigma_vlb * gauss.dev(); PosG[3] *=z2;
			PosG[4] += sigma_vlb * gauss.dev(); PosG[4] *=z2;
			PosG[5] += sigma_vr * gauss.dev();
		}
		else
		if (errorOption == 4)
		{
			if (PosG[2] != 0.0)
				PosG[2] = Star.mag(0)-3.26128+9.41984*(Star.mag(1)-Star.mag(2));
			if (PosG[2] < 0.0)
					PosG[2] = 0.1;

			PosG[3] += sigma_vlb * gauss.dev();
			PosG[4] += sigma_vlb * gauss.dev();
			PosG[5] += sigma_vr * gauss.dev();
		}
		else
		if (errorOption == 5)
		{
			if (PosG[2] != 0.0)
			{
				double z=-1.0;
				while(z < 0.0)
					z= PosG[2]*(1 + sigma_r * gauss.dev());
				PosG[2] =z;
			}
			PosG[3] += sigma_vlb * gauss.dev();
			PosG[4] += sigma_vlb * gauss.dev();
			PosG[5] += sigma_vr * gauss.dev();
		}


		Cot::lbr_to_xyz(PosG, &Star.pos(0));

		Star.feh() += sigma_fe * gauss.dev();
		Star.alpha() += sigma_al * gauss.dev();

	}
}

//
//// 3 kpc 3 km.s  to  rad,rad kpc 3 km/s
//// longitude l --> -Pi to +Pi ,  latitude b --> -Pi/2  to Pi/2
//void SurveyError::xyz_to_lbr(double *PosC,double *PosS)
//{
//    PosS[2]=sqrt(PosC[0]*PosC[0]+PosC[1]*PosC[1]+PosC[2]*PosC[2]);
//    PosS[0]=atan2(PosC[1],PosC[0]);
//    if(PosS[2]>0.0)
//	PosS[1]=asin(PosC[2]/PosS[2]);
//    else
//	PosS[1]=0.0;
//
//    double sin_b,cos_b,sin_l,cos_l;
//    sin_l=sin(PosS[0]);    cos_l=cos(PosS[0]); sin_b=sin(PosS[1]);    cos_b=cos(PosS[1]);
//    TM[0][0]=-sin_l;    TM[0][1]=-sin_b*cos_l;  TM[0][2]= cos_b*cos_l;
//    TM[1][0]= cos_l;    TM[1][1]=-sin_b*sin_l;  TM[1][2]= cos_b*sin_l;
//    TM[2][0]=   0.0;    TM[2][1]= cos_b;        TM[2][2]= sin_b;
//    TM.Tmult((PosC+3),(PosS+3));
//}
//
//
////  rad,rad  kpc 3 km/s to 3 kpc 3 km.s
//// longitude l --> -Pi to +Pi ,  lat b --> -Pi/2  to Pi/2
//void SurveyError::lbr_to_xyz(double *PosS,double *PosC)
//{
//    PosC[0]=PosS[2]*cos(PosS[1])*cos(PosS[0]);
//    PosC[1]=PosS[2]*cos(PosS[1])*sin(PosS[0]);
//    PosC[2]=PosS[2]*sin(PosS[1]);
//    TM.mult((PosS+3),(PosC+3));
//}


SurveyError::~SurveyError()
{
	// TODO Auto-generated destructor stub
}
