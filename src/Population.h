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

#ifndef POPULATION_H_
#define POPULATION_H_

#include "Sampler.h"
#include "StarParticle.h"
#include "GInterpolator.h"




class Population
{
public:
	Population();
	virtual ~Population();
//	virtual ~Population();
	int optionE;
	int ID;
	double a2,Rsun,mass_split,l_split,te[3],PosB[3],PosC[3];
	double rho0,d0,hr1,hr2,eps,rho_fac;
	double age,dage,feh,dfeh;//,dfehdr;
	double sigma_v[3];//,dlnsigmau2dr;
	double vconst;
	Sampler* imfP;
	Interp* vcircP;
	void setParams(int i);
	void setVel(StarParticle &Star);
	void print()
	{
		cout<<left<<setw(12)<<"ID="<<setw(12)<<ID<<setw(12)<<"rho(Rsun)="<<setw(12)<<density(PosC,age)<<setw(12)<<" rho0="<<setw(12)<<rho0<<endl;
		cout<<left<<setw(12)<<"e="<<setw(12)<<eps<<setw(12)<<"Age="<<setw(12)<<age<<setw(12)<<"[FeH]="<<setw(12)<<feh<<endl;
		cout<<left<<setw(12)<<"s_U="<<setw(12)<<sigma_v[0]<<setw(12)<<"s_V="<<setw(12)<<sigma_v[1]<<endl;
//		cout<<left<<setw(12)<<"s_W="<<setw(12)<<sigma_v[2]<<setw(24)<<"Asymmetric Drift km/s="<<asymmetricDrift(PosC1,Rsun,226.0,age,sigma_rr)<<endl;
	}
	inline double pot_eff12(double Rl,double R, double L)
	{
		return log(Rl/R)*vconst*vconst+0.5*L*L*(1/(Rl*Rl)-1/(R*R));
	}
	double vphi_dist(double age1,double R)
	{
		double temp1=0,temp2=1,vphi=0,sigma1,Rl;
		double temp_max=2*surface_density(R)/sigmar_func(age1,R);
		while(temp2>temp1)
		{
			vphi=grandomu()*400.0;
			Rl=vphi*R/vconst;
			sigma1=sigmar_func(age1,Rl);
			temp1=surface_density(Rl)*exp(pot_eff12(Rl,R,vphi*R)/(sigma1*sigma1))/(sigma1);
			temp2=grandomu()*temp_max;
			certify(temp1<=temp_max,"increase max value in vphi_dist");
		}
		return vphi-vconst;
	}
	inline double dlnrho_dlnR(const double* Pos,double age1)
	{
		double d1=0.0,d2=0.0,dr;
		dr=1e-3;
		te[0]=Pos[0]*(1+dr);		te[1]=Pos[1]*(1+dr);		te[2]=Pos[2];
		d2=density(te,age1);
		te[0]=Pos[0]*(1-dr);		te[1]=Pos[1]*(1-dr);		te[2]=Pos[2];
		d1=density(te,age1);
		return (d2-d1)/(dr*(d1+d2));
	}
	inline double radiusC(const double* Pos)
	{
		return sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]);
	}
	inline double radius(const double* Pos)
	{
		return sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]+Pos[2]*Pos[2]);
	}
	inline double radius2(const double* Pos)
	{
		return (Pos[0]*Pos[0]+Pos[1]*Pos[1]+Pos[2]*Pos[2]);
	}
	virtual double eps_func(double age1)
	{
		return eps;
	}
	virtual double surface_density(double R)
	{
		certify(ID<=7,"popID should be less than 7 to use surface density");
		return 1.0;
	}
	virtual double sigmar_func(double age1,double rc)
	{
		return sigma_v[0];
	}
	virtual inline double asymmetricDrift(const double* Pos,double v_c,double age1,double sigma_rr)
	{
		return v_c;
	}
	virtual void computeMinMaxDensity(const double* Pos,double l,double lage,double& rho_min,double& rho_max)
	{
//		cout<<"Population computeMinMaxDensity"<<endl;
		setPosB(Pos);
		double xmind[3],xmaxd[3];
		for(int i=0;i<3;++i)
		{
			if(PosB[i]>0)
			{
				xmaxd[i]=(PosB[i]-l);
				xmind[i]=(PosB[i]+l);
			}
			else
			{
				xmaxd[i]=(PosB[i]+l);
				xmind[i]=(PosB[i]-l);
			}
			if (fabs(PosB[i]) <= l)
//			if((PosB[i]-l)*(PosB[i]+l)<=0)
			{
				xmaxd[i]=0.0;
			}
		}
		rho_min=density(xmind,0.0);
		rho_max=density(xmaxd,0.0);
	}
	virtual double density(const double* Pos,double age1)=0;
	virtual void setPosB(const double* Pos)
	{
		PosB[0]=Pos[0];
		PosB[1]=Pos[1];
		PosB[2]=Pos[2];
	}
	virtual double AMR(const double* Pos,double age1)
	{
		return feh;
	}
	virtual void setAgeFehMass(StarParticle &Star)
	{
		if(optionE==0)
			Star.age()=(age+(grandomu()-0.5)*2*dage);
		Star.feh()=AMR(&Star.pos(0),Star.age())+grandomn()*dfeh*1.0;
		Star.smass() = imfP->rand();
		Star.age()=log10(Star.age());
	}
};


class ThinDisk:public Population
{
public:
	double maxloc;
	vector<double> agevec,fehvec,dfehvec,dagevec,epsvec,sigmauvec;
	double th_max,g_warp,R_warp,R_flare,g_flare,k_flare,sfr;
	double th,R,qparam,bparam;
	GInterpolator surface_density_data;
	ThinDisk(int warpFlareOn1,int optionE1)
	{
		qparam=0.33;
		bparam=0.33;
		optionE=optionE1;
		sfr=2.37;
		if(fabs(Rsun-8.5)<0.01)
			sfr=2.85;

		th_max=3.14159*0.5;
		g_warp=0.18;  		g_flare=5.4e-4; //corrected from 0.54e-4 on April21 2010 from 5.4e-3 Jul29-2009
		R_warp=8.4*Rsun/8.5;         R_flare=9.5*Rsun/8.5;
//		R_flare=1e6;
//		R_warp=1e6;
		{
			agevec.resize(7);
			dagevec.resize(7);
			fehvec.resize(7);
			dfehvec.resize(7);
			epsvec.resize(7);
			double agev[]  ={0.075e9, 0.575e9 ,1.5e9, 2.5e9, 4e9,   6e9,   8.5e9};
			double fehv[]  ={0.01,    0.03,    0.03,  0.01, -0.07, -0.14, -0.37};
			double dfehv[] ={0.12,    0.12,    0.10,  0.11,  0.18,  0.17,  0.2};
			double dagev[] ={0.075e9, 0.425e9, 0.5e9, 0.5e9, 1e9,   1e9,   1.5e9};
			double epsv[] ={0.014,    0.0268,  0.0375,0.0551,0.0696,0.0785,0.0791};
			for(size_t i=0;i<agevec.size();++i)
			{
				agevec[i]=agev[i];
				dagevec[i]=dagev[i];
				fehvec[i]=fehv[i];
				dfehvec[i]=dfehv[i];
				epsvec[i]=epsv[i];
			}
		}
		sigma_v[0]=50.0;
		sigma_v[1]=50.0/1.548;
		sigma_v[2]=50.0/2.38;

		if(warpFlareOn1==0)
		{
			R_warp=1e6;
			R_flare=1e6;
		}
	}
	~ThinDisk();
	inline double eps_func(double age1)
	{
		if(optionE==1)
			return min(0.079,0.104*pow((age1+0.1e9)/10.1e9,0.5));
		else
			return eps;
	}
	double sigmar_func(double age1,double rc)
	{
//		double temp=min(43.1,50.0*pow((age1+0.1e9)/10.1e9,bparam));
//		temp*=pow(density_simple(rc,age1)/density_simple(Rsun,age1),qparam);
//		temp*=pow(surface_density(rc)/surface_density(Rsun),qparam);
//		return temp;
		return min(43.1,50.0*pow((age1+0.1e9)/10.1e9,bparam))*exp((Rsun-rc)*qparam/2.5);
	}
	void setPosB(const double* Pos)
	{
		PosB[0]=Pos[0];
		PosB[1]=Pos[1];
		PosB[2]=Pos[2];
		R=radiusC(Pos);
		th=atan2(Pos[1],Pos[0]);
		if(th<0)
			th=2*3.14159+th;
		if(R>R_warp)
		{
			PosB[2]=Pos[2]-g_warp*(R-R_warp)*cos(th-th_max);
		}
		k_flare=1.0;
		if(R>R_flare)
			k_flare=1+g_flare*(R-R_flare);
	}
	void computeMinMaxDensity(const double* Pos,double l,double lage,double& rho_min,double& rho_max)
	{
//		cout<<"Thin Disk computeMinMaxDensity"<<endl;
		double a_min, a_max, a_nor;
		setPosB(Pos);   //R , th calculated
		//-----------k_flare1(small) k_flare2(large)----------------------------------------
		double k_flare1,k_flare2;
		double eps1,eps2;
		eps1=eps_func(Pos[3]-lage);
		eps2=eps_func(Pos[3]+lage);

		if ((R-l*1.41421) > R_flare)
			k_flare1=1+g_flare*(R-R_flare-l*1.41421);
		else
			k_flare1=1.0;

		if ((R+l*1.41421) > R_flare)
			k_flare2=1+g_flare*(R-R_flare+l*1.41421);
		else
			k_flare2=1.0;

		//---------Warp calculation------------------------------------------
		if ((R+1.41421*l) > R_warp)
		{
			//---------gw and dgw------------------------------------------
			double dth,temp,temp1,temp2;
			temp=1.73205*l/R;
			if(temp>1.0)
				dth=2*3.14159;
			else
				dth = asin(temp);

			temp=(th-th_max);
			if((th-th_max)>3.14159)
				temp=-(3.1415982-temp);
			temp=fabs(temp);
			assert((temp<=3.14159));

			if((temp-dth) < 0)
				temp1=1.0;
			else
				temp1=cos(temp-dth);

			if((temp+dth) > 3.14159)
				temp2=-1.0;
			else
				temp2=cos(temp+dth);
			temp=cos(temp);
			temp2=g_warp*(max((temp1-temp),(temp-temp2)));
			assert(temp2>=0);
			//-----------total dz->temp----------------------------------------
			temp1=g_warp*fabs(cos(th-th_max));
			temp=1.41421*l+temp2*fabs(R-R_warp)+(temp1+temp2)*1.41421*l;
			//-----------a_min a_max----------------------------------------
			temp1=(fabs(PosB[2])+temp)/(eps1*k_flare1);
			a_max=sqrt((R+1.41421*l)*(R+1.41421*l)+temp1*temp1);

			if(R>1.41421*l)
				a_min=(R-1.41421*l)*(R-1.41421*l);
			else
				a_min=0.0;
			temp1=max(fabs(PosB[2])-temp,0.0)/(eps2*k_flare2);
			if(temp>0)
				a_min+=temp1*temp1;
			a_min=sqrt(a_min);
		}
		//------------- no warp-------------------------------------------------------
		else
		{
			double xmind[3], xmaxd[3];
			for (int i = 0; i < 3; ++i)
			{
				if (PosB[i] > 0)
				{
					xmaxd[i] = (PosB[i] - l);
					xmind[i] = (PosB[i] + l);
				}
				else
				{
					xmaxd[i] = (PosB[i] + l);
					xmind[i] = (PosB[i] - l);
				}
//				if (fabs(PosB[i]) <= l)
				if ((PosB[i] - l) * (PosB[i] + l) <= 0)
				{
					xmaxd[i] = 0.0;
				}
			}

			a_min = sqrt(xmaxd[0] * xmaxd[0] + xmaxd[1] * xmaxd[1] + (xmaxd[2]
					* xmaxd[2]) / (eps2 * eps2*k_flare2*k_flare2));
			a_max = sqrt(xmind[0] * xmind[0] + xmind[1] * xmind[1] + (xmind[2]
					* xmind[2]) / (eps1 * eps1*k_flare1*k_flare1));
		}

		a_nor = sqrt(PosB[0] * PosB[0] + PosB[1] * PosB[1] + (PosB[2]
				* PosB[2]) / ((eps1+eps2) * (eps1+eps2)*0.25*k_flare*k_flare));

		if (a_nor > maxloc)
		{
			rho_min = density_simple(a_max,Pos[3]+lage)/(k_flare2);
			rho_max = density_simple(a_min,Pos[3]-lage)/(k_flare1);
		}
		else
		{
			rho_min = density_simple(a_min,Pos[3]+lage)/(k_flare2);
			rho_max = density_simple(a_max,Pos[3]-lage)/(k_flare1);
		}
		if ((a_max >= maxloc) && (a_min <= maxloc))
			rho_max = density_simple(maxloc,Pos[3]-lage)/k_flare1;
	}
//	void computeMinMaxDensity(const double* Pos,double l,double& rho_min,double& rho_max)
//	{
////		cout<<"Thin Disk computeMinMaxDensity"<<endl;
//		double a_min, a_max, a_nor;
//		setPosB(Pos);   //R , th calculated
//		//-----------k_flare1(small) k_flare2(large)----------------------------------------
//		{
//			double xmind[3], xmaxd[3];
//			for (int i = 0; i < 3; ++i)
//			{
//				if (PosB[i] > 0)
//				{
//					xmaxd[i] = (PosB[i] - l);
//					xmind[i] = (PosB[i] + l);
//				}
//				else
//				{
//					xmaxd[i] = (PosB[i] + l);
//					xmind[i] = (PosB[i] - l);
//				}
//				if (fabs(PosB[i]) <= l)
////				if ((PosB[i] - l) * (PosB[i] + l) <= 0)
//				{
//					xmaxd[i] = 0.0;
//				}
//			}
//
//			a_min = sqrt(xmaxd[0] * xmaxd[0] + xmaxd[1] * xmaxd[1] + (xmaxd[2]
//					* xmaxd[2]) / (eps * eps));
//			a_max = sqrt(xmind[0] * xmind[0] + xmind[1] * xmind[1] + (xmind[2]
//					* xmind[2]) / (eps * eps));
//		}
//
//		a_nor = sqrt(PosB[0] * PosB[0] + PosB[1] * PosB[1] + (PosB[2]
//				* PosB[2]) / (eps * eps*k_flare*k_flare));
//
//		if (a_nor > maxloc)
//		{
//			rho_min = density_simple(a_max);
//			rho_max = density_simple(a_min);
//		}
//		else
//		{
//			rho_min = density_simple(a_min);
//			rho_max = density_simple(a_max);
//		}
//		if ((a_max >= maxloc) && (a_min <= maxloc))
//			rho_max = density_simple(maxloc);
//	}
	double density(const double* Pos,double age1)
	{
		return 0.0;
	}
	double AMR(const double* Pos,double age1)
	{
		return interpolate(agevec,fehvec,age1)-0.07*(sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1])-Rsun);
	}
	double asymmetricDrift(const double* Pos,double v_c,double age1,double sigma_rr)
	{
		double temp=0;//,v_c,r_c;
		double r_c=radiusC(Pos);
		temp=(1-sigma_v[1]*sigma_v[1]/(sigma_v[0]*sigma_v[0]))+(1-sigma_v[2]*sigma_v[2]/(sigma_v[0]*sigma_v[0])) ;
		temp+=(dlnrho_dlnR(Pos,age1)-2*qparam*r_c/2.5);
//		temp+=(dlnrho_dlnR(Pos,age1)+2*qparam*dlnrhoSurf_dlnR(r_c));
//		temp+=(1+2*qparam)*(dlnrho_dlnR(Pos,age1));
		temp*=(sigma_rr*sigma_rr*0.5/v_c);
		temp=-temp; //change sign
		return temp;
	}
	double surface_density(double r)
	{
		return rho0*surface_density_data.interpol(&r);
	}
private:
//	inline double dlnrhoSurf_dlnR(double rc)
//	{
//		double d1=0.0,d2=0.0,dr;
//		dr=1e-3;
//		d1=surface_density(rc*(1-dr));
//		d2=surface_density(rc*(1+dr));
//		return (d2-d1)/(dr*(d1+d2));
//	}
	virtual double density_simple(double a,double age1)=0;
};



class ThinDisk0:public ThinDisk
{
public:
//	double maxloc;
//	double th_max,g_warp,R_warp,R_flare,g_flare,k_flare;
//	double th,R,theta;
	ThinDisk0(int warpFlareOn1,int optionE1,const string& inputDir,Interp* vcircP1):ThinDisk(warpFlareOn1,optionE1)
	{
		vcircP=vcircP1;
		ID=0;
		age=agevec[ID];
		dage=dagevec[ID];
		feh=fehvec[ID];
		dfeh=dfehvec[ID];
		eps=epsvec[ID];
		setParams(ID);
//		th_max=90.0*3.14159/180.0;		g_warp=0.18; R_warp=1000;
//		R_flare=1000.0;		g_flare=5.4e-4;
		maxloc=3.79;
		hr1=5.0;
		hr2=3.0;
		rho0=sfr/(5.5683282*(hr1*hr1*hr1-hr2*hr2*hr2));
		if(optionE1==0)
			rho0*=(2*dage);
//		rho0=sfr*2*dage/(5.5683282*(hr1*hr1*hr1-hr2*hr2*hr2));
//		dage=0.0;
		hr1*=hr1;
		hr2*=hr2;
		surface_density_data.readFromFile(inputDir+"Model/surface_density.ebf","/Surface_density1");
//		d0=density(PosC)/rho0;
//		rho0=rho0/d0;
//		rho0=sfr*2*dage/(5.5683282*(hr1*hr1*hr1-hr2*hr2*hr2));
//		cout<<"rho0="<<rho0/eps_func(age)<<endl;
//		dage=0.0;

	}
	~ThinDisk0();
	double density_simple(double a,double age1)
	{
		return (rho0/eps_func(age1))*(exp(-a*a/hr1)-exp(-a*a/hr2));
	}
	double density(const double* Pos,double age1)
	{
		setPosB(Pos);
		double eps1=eps_func(age1);
		a2=PosB[0]*PosB[0]+PosB[1]*PosB[1]+(PosB[2]*PosB[2])/(eps1*eps1*k_flare*k_flare);
//		if((a2*c1)<0.75)
//		{
//			double temp,temp1;
//			temp=0.0;
//			temp1=1.0;
//			for(int i=1;i<5;++i)
//			{
//				temp1*=(a2*c1/i);
//				temp+=temp1;
//			}
//			return (rho0)*temp;
//		}
		return (rho0/eps1)*(exp(-a2/hr1)-exp(-a2/hr2))/(k_flare);
	}
};


class ThinDisk1:public ThinDisk
{
public:
//	double maxloc;
//	double th_max,g_warp,R_warp,R_flare,g_flare,k_flare;
//	double th,R,theta;
	ThinDisk1(int i,int warpFlareOn1,int optionE1,const string& inputDir,Interp* vcircP1):ThinDisk(warpFlareOn1,optionE1)
	{
		vcircP=vcircP1;
		ID=i;
		age=agevec[ID];
		dage=dagevec[ID];
		feh=fehvec[ID];
		dfeh=dfehvec[ID];
		eps=epsvec[ID];
		setParams(ID);
//		th_max=90.0*3.14159/180.0;		g_warp=0.18; R_warp=1000;
//		R_flare=1000.0;		g_flare=5.4e-4;
		maxloc=0.229;
		hr1=2.530;
		hr2=1.320;
		rho0=sfr/(23.719602*(hr1*hr1*hr1-hr2*hr2*hr2));
		if(ID==6)
			rho0=rho0/1.25;
		if(optionE1==0)
			rho0*=(2*dage);
//		dage=0.0;
		hr1*=hr1;
		hr2*=hr2;
		surface_density_data.readFromFile(inputDir+"Model/surface_density.ebf","/Surface_density2");

//		d0=density(PosC)/rho0;
//		rho0=rho0/d0;
//		cout<<"rho0="<<rho0<<endl;
	}
	~ThinDisk1();
	double density_simple(double a,double age1)
	{
		double eps1=eps_func(age1);
		return (rho0/eps1)*(exp(-sqrt(0.25+a*a/(hr1)))-exp(-sqrt(0.25+a*a/(hr2))));
	}
	double density(const double* Pos,double age1)
	{
		setPosB(Pos);
		double eps1=eps_func(age1);
		a2=PosB[0]*PosB[0]+PosB[1]*PosB[1]+(PosB[2]*PosB[2])/(eps1*eps1*k_flare*k_flare);
		return (rho0/eps1)*(exp(-sqrt(0.25+a2/hr1))-exp(-sqrt(0.25+a2/hr2)))/(k_flare);
	}
};



class ExpDisk:public Population
{
public:
	ExpDisk()
	{
		optionE=0;
		hr1=3.5;
		hr2=0.2;
		eps=0.0268;
		d0=0.0881;
		rho0=7.9e6;
	}
	~ExpDisk();
	double density(const double* Pos,double age1)
	{
		a2=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]);
		return (rho0/d0)*exp(-(a2/hr1+fabs(Pos[2])/hr2));
	}
};

class ThickDisk:public Population
{
public:
	double xl,c1,c2;
	double th_max,g_warp,R_warp,R_flare,g_flare,k_flare;
	double th,R,qparam,bparam;
	ThickDisk(int warpFlareOn1,Interp* vcircP1 )
	{
		optionE=0;
		qparam=0.33;
		bparam=0.33;

		th_max=3.14159*0.5;
		g_warp=0.18;  		g_flare=5.4e-4; //corrected from 0.54e-4 on April21 2010
		R_warp=8.51*Rsun/8.5;         R_flare=9.5*Rsun/8.5;

		if(warpFlareOn1==0)
		{
			R_flare=1e6;
			R_warp=1e6;
		}

		age=11e9;
		dage=0.0;
		feh=-0.78;
		dfeh=0.3;
//		dfeh=0.0;

		eps=0.0;
		sigma_v[0]=67.0;
		sigma_v[1]=51.0;
		sigma_v[2]=42.0;

		vcircP=vcircP1;
		ID=7;
		setParams(ID);
		hr1=2.5;
		hr2=0.8;
		xl=0.4;
		c1=1/((hr2*xl)*(2+xl/hr2));
		c2=exp(xl/hr2)/(1+0.5*xl/hr2);
		rho0=1.34e6*1157.68;
		d0=density(PosC,age)/rho0;
		rho0=rho0/d0;

	}
	~ThickDisk();
	void setPosB(const double* Pos)
	{
		PosB[0]=Pos[0];
		PosB[1]=Pos[1];
		PosB[2]=Pos[2];
		R=radiusC(Pos);
		th=atan2(Pos[1],Pos[0]);
		if(th<0)
			th=2*3.14159+th;
		if(R>R_warp)
		{
			PosB[2]=Pos[2]-g_warp*(R-R_warp)*cos(th-th_max);
		}
		k_flare=1.0;
		if(R>R_flare)
			k_flare=1+g_flare*(R-R_flare);
	}
	double density(const double* Pos,double age1)
	{
		setPosB(Pos);
		a2=sqrt(PosB[0]*PosB[0]+PosB[1]*PosB[1]);
//		if(fabs(PosB[2])<xl*k_flare)
//			return rho0*exp(-(a2-Rsun)/hr1)*(1-c1*PosB[2]*PosB[2]/(k_flare*k_flare))/k_flare;
//		return rho0*exp(-(a2-Rsun)/hr1)*c2*exp(-fabs(PosB[2])/(hr2*k_flare))/k_flare;

		double temp=(2*hr2+xl)*(3*k_flare*hr2*xl+3*k_flare*k_flare*hr2*hr2+xl*xl)/(
				(2*hr2*k_flare+xl)*(3*hr2*xl+3*hr2*hr2+xl*xl));
		if(fabs(Pos[2])<xl)
			return rho0*exp(-(a2-Rsun)/hr1)*(1-PosB[2]*PosB[2]/(xl*(2*hr2*k_flare+xl)))/temp;
		return rho0*exp(-(a2-Rsun)/hr1)*(exp(xl/(hr2*k_flare))/(1+0.5*xl/(hr2*k_flare)))*exp(-fabs(PosB[2])/(hr2*k_flare))/temp;


	}
	void computeMinMaxDensity(const double* Pos,double l,double lage,double& rho_min,double& rho_max)
	{
//		cout<<"Thin Disk computeMinMaxDensity"<<endl;
		setPosB(Pos);   //R , th calculated

		double xmind[3], xmaxd[3];
		for (int i = 0; i < 3; ++i)
		{
			if (PosB[i] > 0)
			{
				xmaxd[i] = (PosB[i] - l);
				xmind[i] = (PosB[i] + l);
			}
			else
			{
				xmaxd[i] = (PosB[i] + l);
				xmind[i] = (PosB[i] - l);
			}
//				if (fabs(PosB[i]) <= l)
			if ((PosB[i] - l) * (PosB[i] + l) <= 0)
			{
				xmaxd[i] = 0.0;
			}
		}

		//-----------k_flare1(small) k_flare2(large)----------------------------------------
		double k_flare1,k_flare2;

		if ((R-l*1.41421) > R_flare)
			k_flare1=1+g_flare*(R-R_flare-l*1.41421);
		else
			k_flare1=1.0;

		if ((R+l*1.41421) > R_flare)
			k_flare2=1+g_flare*(R-R_flare+l*1.41421);
		else
			k_flare2=1.0;

		//---------Warp calculation------------------------------------------
		if ((R+1.41421*l) > R_warp)
		{
			//---------gw and dgw------------------------------------------
			double dth,temp,temp1,temp2;
			temp=1.73205*l/R;
			if(temp>1.0)
				dth=2*3.14159;
			else
				dth = asin(temp);

			temp=(th-th_max);
			if((th-th_max)>3.14159)
				temp=-(3.1415982-temp);
			temp=fabs(temp);
			assert((temp<=3.14159));

			if((temp-dth) < 0)
				temp1=1.0;
			else
				temp1=cos(temp-dth);

			if((temp+dth) > 3.14159)
				temp2=-1.0;
			else
				temp2=cos(temp+dth);
			temp=cos(temp);
			temp2=g_warp*(max((temp1-temp),(temp-temp2)));
			assert(temp2>=0);
			//-----------total dz->temp----------------------------------------
			temp1=g_warp*fabs(cos(th-th_max));
			temp=1.41421*l+temp2*fabs(R-R_warp)+(temp1+temp2)*1.41421*l;
			//-----------a_min a_max----------------------------------------
			xmind[2]=(fabs(PosB[2])+temp);
			xmaxd[2]=max(fabs(PosB[2])-temp,0.0);
		}
		rho_min = min(density_simple(xmind,k_flare2),density_simple(xmind,k_flare1));
		rho_max = max(density_simple(xmaxd,k_flare1),density_simple(xmaxd,k_flare2));
	}
	double sigmar_func(double age1,double rc)
	{
		return sigma_v[0];
//		return sigma_v[0]*exp((Rsun-rc)*qparam/2.5);
	}
	double surface_density(double R)
	{
		return rho0*exp(-R/hr1);
	}
	double asymmetricDrift(const double* Pos,double v_c,double age1,double sigma_rr)
	{
		double temp=0,r_c;//,v_c,r_c;
		r_c=radiusC(Pos);
		temp=(1-sigma_v[1]*sigma_v[1]/(sigma_v[0]*sigma_v[0]))+(1-sigma_v[2]*sigma_v[2]/(sigma_v[0]*sigma_v[0])) ;
		temp-=r_c*(1.0/hr1+(2*qparam)/2.5);
//		temp-=(1+2*qparam)*sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1])/hr1;
		temp*=(sigma_rr*sigma_rr*0.5/v_c);
		temp=-temp; //change sign
		return temp;
	}
private:
	double density_simple(const double* Pos,double k_flarec)
	{
//		a2=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]);
//		if(fabs(Pos[2])<xl*k_flarec)
//			return rho0*exp(-(a2-Rsun)/hr1)*(1-c1*Pos[2]*Pos[2]/(k_flarec*k_flarec))/k_flarec;
//		return rho0*exp(-(a2-Rsun)/hr1)*c2*exp(-fabs(Pos[2])/(hr2*k_flarec))/k_flarec;
		double temp=(2*hr2+xl)*(3*k_flarec*hr2*xl+3*k_flarec*k_flarec*hr2*hr2+xl*xl)/(
				(2*hr2*k_flarec+xl)*(3*hr2*xl+3*hr2*hr2+xl*xl));
		a2=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]);
		if(fabs(Pos[2])<xl)
			return rho0*exp(-(a2-Rsun)/hr1)*(1-Pos[2]*Pos[2]/(xl*(2*hr2*k_flarec+xl)))/temp;
		return rho0*exp(-(a2-Rsun)/hr1)*(exp(xl/(hr2*k_flarec))/(1+0.5*xl/(hr2*k_flarec)))*exp(-fabs(Pos[2])/(hr2*k_flarec))/temp;
	}

};


class Bulge:public Population
{
public:
	double a,b,g,x0,y0,z0,Rc;
	Matrix<double> TM;
	Bulge(Interp* vcircP1):TM(3,3)
	{
		optionE=0;
		age=10e9;
		dage=0.0;
		feh=0.0;
		dfeh=0.4;
		eps=0.0;
//		sigma_v[0]=113.0;
//		sigma_v[1]=115.0;
//		sigma_v[2]=100.0;

		sigma_v[0]=110.0;
		sigma_v[1]=110.0;
		sigma_v[2]=100.0;

		vcircP=vcircP1;
//		double PI=3.141392;
		ID=9;
		setParams(ID);
		rho0=0.255*13.7e9;
		a=78.9 ; b=3.5; g=91.3;

		x0=1.59; y0=0.424; z0=0.424;
		Rc=2.54;
		TM=Cot::RotationMatrix(0,(g)*PI/180.0);
		TM*=Cot::RotationMatrix(1,(-b)*PI/180.0);
		TM*=Cot::RotationMatrix(2,(a-90.0)*PI/180.0);
		//		TM=Matrix<double>(3,3);
		//		Matrix<double> TR1(3,3);
		//		TM.identity();
	}
	~Bulge();
	void setPosB(const double* Pos)
	{
		TM.mult(Pos,PosB);
	}
	void computeMinMaxDensity(const double* Pos,double l,double lage,double& rho_min,double& rho_max)
	{
//		cout<<"Population computeMinMaxDensity"<<endl;
		double l1=l*1.7328;
		setPosB(Pos);
		double xmind[3],xmaxd[3];
		for(int i=0;i<3;++i)
		{
			if(PosB[i]>0)
			{
				xmaxd[i]=(PosB[i]-l1);
				xmind[i]=(PosB[i]+l1);
			}
			else
			{
				xmaxd[i]=(PosB[i]+l1);
				xmind[i]=(PosB[i]-l1);
			}
			if (fabs(PosB[i]) <= l1)
			{
				xmaxd[i]=0.0;
			}
		}
		rho_min=density_simple(xmind,age);
		rho_max=density_simple(xmaxd,age);
	}
	double density(const double* Pos,double age1)
	{
		setPosB(Pos);
		double temp;
		temp=PosB[0]*PosB[0]/(x0*x0)+PosB[1]*PosB[1]/(y0*y0);
		a2=temp*temp;
		temp=PosB[2]*PosB[2]/(z0*z0);
		a2=sqrt(a2+temp*temp);
		temp=sqrt(PosB[0]*PosB[0]+PosB[1]*PosB[1])-Rc;
		if(temp < 0)
			return  rho0*exp(-0.5*a2);
		return  rho0*exp(-0.5*a2)*exp(-temp*temp/0.5);
	}
	double asymmetricDrift(const double* Pos,double v_c,double age1,double sigma_rr)
	{
		return (v_c-radiusC(Pos)*71.62);
//		return 79.0;
	}
private:
	double density_simple(const double* Pos,double age1)
	{
		double temp;
		temp=Pos[0]*Pos[0]/(x0*x0)+Pos[1]*Pos[1]/(y0*y0);
		a2=temp*temp;
		temp=Pos[2]*Pos[2]/(z0*z0);
		a2=sqrt(a2+temp*temp);
		temp=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1])-Rc;
		if(temp < 0)
			return  rho0*exp(-0.5*a2);
		return  rho0*exp(-0.5*a2)*exp(-temp*temp/0.5);
	}
};


class Spheroid:public Population
{
public:
	double ac,rhoc;
	Spheroid(Interp* vcircP1)
	{
		optionE=0;
		age=14e9;
//		age=13.0e9;
		dage=0.0;
		feh=-1.78;
		dfeh=0.5;
		eps=0.76;
//		eps=0.64;
//		sigma_v[0]=131.0;
//		sigma_v[1]=106.0;
//		sigma_v[2]=85.0;
		sigma_v[0]=135.0;
		sigma_v[1]=85.0;
		sigma_v[2]=85.0;

		vcircP=vcircP1;
		ID=8;
		setParams(ID);
		ac=0.5;
		rho0=9.32e3*1407.52*1.1;
//		rho0=9.32e3*1407.52*2.15;
		d0=density(PosC,age)/rho0;
		rho0=rho0/d0;
		rhoc=rho0*pow(ac/Rsun,-2.44);
	}
	~Spheroid();
	double density(const double* Pos,double age1)
	{
		a2=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]+(Pos[2]*Pos[2])/(eps*eps));
		if(a2<ac)
			return  rhoc;
		return  rho0*pow(a2/Rsun,-2.44);
	}
};


#endif /* POPULATION_H_ */
