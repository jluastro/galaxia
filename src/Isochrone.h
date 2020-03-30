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

#ifndef ISOCHRONE_H_
#define ISOCHRONE_H_

#include "StarParticle.h"
class Isochrone
{
public:
	Isochrone();
	 ~Isochrone();
	 void print();
	 void setMon();
	 void setTip();
	 void setmaxV(double maxV1);
	 void interpolate(StarParticle &Star);
	 void interpolateTGM(float smass,vector<double> &x);
	 double age,FeH,alpha;
	void addDwarfs(int nmags,int dwarfOn)
	{
		double temp[]={0.010,0.012,0.015,0.020,0.030,0.040,0.050,0.055,0.060,0.070,0.072,0.075,0.080,0.090,0.100};
		if(dwarfOn==1)
		{
		m.assign(temp,temp+sizeof(temp)/sizeof(temp[0]));
		Mact.assign(temp,temp+sizeof(temp)/sizeof(temp[0]));
		Mag0.resize(sizeof(temp)/sizeof(temp[0]),1e6);
		Mag1.resize(sizeof(temp)/sizeof(temp[0]),1e6);
		Mag2.resize(sizeof(temp)/sizeof(temp[0]),1e6);
		Teff.resize(sizeof(temp)/sizeof(temp[0]),0.0);
		Grav.resize(sizeof(temp)/sizeof(temp[0]),0.0);
		Lum.resize(sizeof(temp)/sizeof(temp[0]),-1e6);
		}

		for(int i=0;i<nmags;++i)
			Mags.push_back(Mag0);
//		cout<<" m size "<<m.size()<<endl;
	}
	vector<double> m;
	vector<double> Mag0;
	vector<double> Mag1;
	vector<double> Mag2;
	vector<double> Teff;
	vector<double> Grav;
	vector<double> Lum;
	vector<double> Mact;
	vector<vector<double> > Mags;
	vector<double> Vmon;
	double m_min,m_max,maxV,mtip;
};

//void setTip1();


#endif /*ISOCHRONE_H_*/
