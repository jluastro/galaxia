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

#include "StellarPopulation.h"

StellarPopulation::StellarPopulation(int i,const double* posC,int warpFlareOn1,Interp* vcircP1,int option,const string &inputDir)
{
	assert((i>=0)&&(i<=9));
	switch (i)
	{
	case 0:
		cpop=new ThinDisk0(warpFlareOn1,1,inputDir,vcircP1);
		cout<<"Thin Disk,ID=0:"<<endl;
		break;
	case 7:
		cpop=new ThickDisk(warpFlareOn1,vcircP1);
		cout<<"ThickDisk:"<<endl;
		break;
	case 8:
		cpop=new Spheroid(vcircP1);
		cout<<"Spheroid:"<<endl;
		break;
	case 9:
		cpop=new Bulge(vcircP1);
		cout<<"Bulge:"<<endl;
		break;
	default:
		cpop=new ThinDisk1(i,warpFlareOn1,1,inputDir,vcircP1);
		cout<<"Thin Disk,ID="<<i<<":"<<endl;
//		cout<<"Thin Disk:"<<i<<endl;
		break;
	}
//	cpop->print();
//	cout<<"rho0"<<" "<<density(PosC)<<" "<<cpop->age<<" "<<cpop->feh<<endl;
//	cout<<"rho0"<<" "<<density(PosC)<<std::endl;

	BHT.initialize(*cpop,posC,warpFlareOn1,option,inputDir);
	Star.popID()=i;
	Star.satID()=0;
	// TODO Auto-generated constructor stub

}

StellarPopulation::~StellarPopulation()
{
	delete cpop;
	// TODO Auto-generated destructor stub
}


void StellarPopulation::spawn(SurveyDesign &sur,IsochroneDB &ic,double fSample)
{
	int ntot1=0;
	//	cout<<"Entering Spawn.. "<<endl;
	ic.setimf(cpop->imfP);
//	double distMod;
	cout<<left<<setw(12)<<"Completed %"<<"<";
	for(size_t i=0;i<BHT.leafs.size();++i)
	{
		BHT.leafs[i]->scale_mass(cpop->rho_fac);
		for(int j=0;j<3;++j)
			Star.pos(j)=BHT.leafs[i]->Pos[j]-sur.posC[j];

		double rc=sqrt(Star.pos(0)*Star.pos(0)+Star.pos(1)*Star.pos(1)+Star.pos(2)*Star.pos(2));
		rc=rc-BHT.leafs[i]->l*1.732;

		if((sur.geo->checkSphere(&Star.pos(0),BHT.leafs[i]->l*1.732)==1)&&(rc < sur.r_max))
		{
			//		double temp=BHT.leafs[i]->mass*1e-6;
		//		int n=int(temp);
		//		if (randomu()<(temp-n)) n++;
			double rc=sqrt(BHT.leafs[i]->Pos[0]*BHT.leafs[i]->Pos[0]+BHT.leafs[i]->Pos[1]*BHT.leafs[i]->Pos[1]);
			double feh_temp=(cpop->feh);
			if(cpop->ID <7)
				feh_temp-=0.07*(rc-cpop->Rsun);
//		ic.setIsochrones(cpop->age,cpop->dage,feh_temp,cpop->dfeh,min(sur.absMag[1],sur.appMag[1]- BHT.leafs[i]->dmod),BHT.leafs[i]->mass*fSample);
		ic.setIsochrones1(cpop->age,cpop->dage,feh_temp,cpop->dfeh,min(sur.absMag[1],sur.appMag[1]- BHT.leafs[i]->dmod),BHT.leafs[i]->mass*fSample,sur.All->starType);

		if((ic.n_tot==1)&&(BHT.leafs[i]->mass<1e-2)&&(BHT.leafs[i]->rho_max/BHT.leafs[i]->rho_min>250.0))
			ic.n_tot=0;

//		ic.n_tot=0;

//		if(i==63781)
//		if(i*1.0/BHT.leafs.size() > 0.3)
//		{
//			cout<<i<<" ntot "<<ic.n_tot<<" defeh "<<cpop->dfeh<<" dage "<<cpop->dage<<" "<<cpop->rho_fac<<endl;
//			BHT.leafs[i]->print();
//		}
		int ntot2=0,ntot3=0;
		for(int j=0;j<ic.n_tot;++j)
		{

//						cout<<i<<" "<<j<<" "<<ic.n_tot<<" "<<cpop->dfeh<<" "<<cpop->dage<<endl;

//			BHT.leafs[i]->generatePos(Star.Pos,*cpop);
			if(cpop->optionE==1)
				BHT.leafs[i]->generatePos(Star,*cpop);

//			ic.generateStar(Star);
			cpop->setAgeFehMass(Star);
			ic.interpolateStar(Star);
			Star.partID() = 0;
			if(Star.mag(0) <1e4)
			{
				if(cpop->optionE==0)
					BHT.leafs[i]->generatePos(Star,*cpop);
				cpop->setVel(Star);
			for(int k=0;k<6;++k)
				Star.pos(k)-=sur.posC[k];

			Star.AbsToAppMag();
			ntot3+=sur.push_back(Star);
			ntot1++;
			ntot2++;
			}

		}
//		if(((ntot2*1.0/ntot3)>10)&&(ntot3>0))
//		{
//			cout<<ntot2<<" "<<ntot3<<endl;
//			cout<<i<<" "<<ic.n_tot/(BHT.leafs[i]->mass*fSample)<<" "<<1.0/cpop->imfP->meanx<<" "<<BHT.leafs[i]->mass<<" "<<ic.n_tot<<endl;
//			cout<<"m_min m_max "<<ic.box_mmin<<" "<<ic.box_mmax<<" "<<(cpop->imfP->getFac(ic.box_mmin,ic.box_mmax)/cpop->imfP->meanx)<<endl;
//		}

		}
		if((i%(BHT.leafs.size()/10))==0)
			cout<<i*100/BHT.leafs.size()<<".."; cout<<flush;
	}

	cout<<">"<<endl;
//	cout<<"ntot1 "<<ntot1<<endl;
}
