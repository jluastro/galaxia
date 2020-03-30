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

#ifndef ISOCHRONEDB_H_
#define ISOCHRONEDB_H_

#include "IsochroneBase.h"
#include "Sampler.h"


class Erf {
private:
	int ncof;
	vector<double> cof;
public:
	Erf():ncof(28)
	{
		double cof1[28] = {-1.3026537197817094, 6.4196979235649026e-1,
			1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
			3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
			-1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
			6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
			9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
			-1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};
		cof.resize(28,0);
		for(size_t i=0;i<cof.size();++i)
			cof[i]=cof1[i];

	}
	inline double erf(double x) {
		if (x >=0.) return 1.0 - erfccheb(x);
		else return erfccheb(-x) - 1.0;
	}

	inline double erfc(double x) {
		if (x >= 0.) return erfccheb(x);
		else return 2.0 - erfccheb(-x);
	}

	double erfccheb(double z){
		int j;
		double t,ty,tmp,d=0.,dd=0.;
		if (z < 0.) throw("erfccheb requires nonnegative argument");
		t = 2./(2.+z);
		ty = 4.*t - 2.;
		for (j=ncof-1;j>0;j--) {
			tmp = d;
			d = ty*d - dd + cof[j];
			dd = tmp;
		}
		return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
	}

	double inverfc(double p) {
		double x,err,t,pp;
		if (p >= 2.0) return -100.;
		if (p <= 0.0) return 100.;
		pp = (p < 1.0)? p : 2. - p;
		t = sqrt(-2.*log(pp/2.));
		x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
		for (int j=0;j<2;j++) {
			err = erfc(x) - pp;
			x += err/(1.12837916709551257*exp(-(x*x))-x*err);
		}
		return (p < 1.0? x : -x);
	}

	inline double inverf(double p) {return inverfc(1.-p);}
	double erfcc(const double x)
	{
		double t,z=fabs(x),ans;
		t=2./(2.+z);
		ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
			t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
			t*(-0.82215223+t*0.17087277)))))))));
		return (x >= 0.0 ? ans : 2.0-ans);
	}


};


class Normaldist : public Erf
{
public:
	double mu, sig;
	Normaldist(double mmu = 0., double ssig = 1.) : mu(mmu), sig(ssig) {
		if (sig <= 0.) throw("bad sig in Normaldist");
	}
	double p(double x) {
		return (0.398942280401432678/sig)*exp(-0.5*((x-mu)*(x-mu)/(sig*sig)));
	}
	double cdf(double x) {
		return 0.5*erfc(-0.707106781186547524*(x-mu)/sig);
	}
	double invcdf(double p) {
		if (p <= 0. || p >= 1.) throw("bad p in Normaldist");
		return -1.41421356237309505*sig*inverfc(2.*p)+mu;
	}
};


class IsochroneDB: public IsochroneBase
{
public:
	IsochroneDB(const string& inputDir,const string &dirname,const string &photoSys,const string &magcolorNames,int extraFieldsOn1,Sampler* imfP1=NULL):IsochroneBase(inputDir,dirname,photoSys,extraFieldsOn1,magcolorNames)
	{
		imfP=imfP1;
		maxV=30;
		delta_maxV=0.01;
		W_Age.resize(Age.size(),0);
		// corrected 6/5/2010  W_FeH.resize(Age.size(),0);
		W_FeH.resize(FeH.size(),0);
		dAge=Age[1]-Age[0];
		dFeH=(FeH.back()-FeH.front())/(FeH.size()-1);
		last_age=0.0;
		last_dage=0.0;
	}
	virtual ~IsochroneDB();
	void setIsochrones(double Age2,double dAge2,double FeH2,double dFeH2,double maxV1,double mass);
	void setIsochrones1(double Age2,double dAge2,double FeH2,double dFeH2,double maxV1,double mass,int starType);
	void generateStar(StarParticle &Star);
	void interpolateStar(StarParticle &Star);
	void interpolateTGM(float age1,float feh1,float smass1,vector<double> &x);
	void min_max_m(vector<double> &age_a, vector<double> &FeH_a,
			vector<double> &Alpha_a, vector<double> &m1_a, vector<double> &m2_a,
			double maxV, int OptRR);
	void setimf(Sampler* imfP1)
	{
		imfP=imfP1;
	}
	int n_tot;
	double box_mmin,box_mmax;
private:
	void calculateAgeProb(double Age2,double dAge2);
	void calculateFeHProb(double FeH2,double dFeH2);
	Sampler* imfP;
	double dAge,dFeH,dAlpha,last_age,last_dage;
	vector<double> W_Age,W_FeH,ic_p;
	double maxV,delta_maxV,p_tot;
	Normaldist gaussian;
};

#endif /* ISOCHRONEDB_H_ */
