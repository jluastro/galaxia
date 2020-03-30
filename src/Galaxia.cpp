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
#include"Parameters.h"
#include "Timer.h"
#include "ebfvector.hpp"
//#include "GInterpolator.h"
#include "Satellite.h"
//#include <stdint.h>

double Chabrier_exp_imf1(double x)
{
	return 22.8978*(exp(-pow((716.4/x),0.25)))*pow(x,-3.3);
}

double Kroupa_imf(double x)
{
	if(x >0.5)
		return pow(x,-2.3);
	else if (x <0.08)
		return pow(x,-0.3);
	else
		return pow(x,-1.3);
}



void append_ext(const string& filename,const string dir)
{
//	double PI=3.14159;
//	double DTOR=0.0174533;
	ebf::EbfDataInfo dinfo;
	double RADEG=57.2958;
	GInterpolator gi;
	gi.readFromFile(dir+"Extinction/ExMap3d_1024.ebf","/ExMap3d");
	GInterpolator ge_solar;
	ge_solar.readFromFile(dir+"Extinction/ExMap3d_1024.ebf","/ExMap2d");
	GInterpolator ge_schlegel;
	ge_schlegel.readFromFile(dir+"Extinction/Schlegel_4096.ebf","/ExMap2d");
	double PosC[3],PosS[3];
//	fbvector<float> x(filename+".ebf","r","/Pos3");
	ebf::EbfVector<float> x(filename+".ebf","/px");
	ebf::EbfVector<float> y(filename+".ebf","/py");
	ebf::EbfVector<float> z(filename+".ebf","/pz");
	vector<float> exbv1(x.dim(0),0.0);
	vector<float> exbv2(x.dim(0),0.0);
	vector<float> exbv3(x.dim(0),0.0);
	vector<float> exbv4(x.dim(0),0.0);
	vector<float> exbv5(x.dim(0),0.0);
//	certify(x.dim(1)==3,"Check map dim");
	for(int64_t i=0;i<x.dim(0);++i)
	{
//		cout<<i<<" "<<x.dim(0)<<" "<<x.size()<<" "<<i*x.dim(0)+2<<endl;
//		PosC[0]=x[i*x.dim(1)+0];	PosC[1]=x[i*x.dim(1)+1];	PosC[2]=x[i*x.dim(1)+2];
		PosC[0]=x[i];	PosC[1]=y[i];	PosC[2]=z[i];
		Cot::xyz_to_lbr(PosC,PosS,3);
		PosS[0]*=RADEG;		PosS[1]*=RADEG;
		if(PosS[0]<0)
			PosS[0]=360.0+PosS[0];
		PosS[2]=log10(PosS[2]);
		exbv1[i]=gi.interpol(PosS);
		exbv2[i]=exbv1[i];
		exbv1[i]*=ge_schlegel.interpol(PosS);
		exbv2[i]*=ge_solar.interpol(PosS);
		exbv3[i]=ge_schlegel.interpol(PosS);
		exbv4[i]=PosS[0];
		exbv5[i]=PosS[1];
	}

	if(x.dim(0) >0)
	{
		std::string filename1=filename+".ebf";
		if(ebf::ContainsKey(filename1,"/ExBV_Schlegel",dinfo)==0)
		{
			ebf::Write(filename1,"/ExBV_Schlegel",&exbv1[0],"a","",exbv1.size());
			ebf::Write(filename1,"/ExBV_Solar",&exbv2[0],"a","",exbv2.size());
			ebf::Write(filename1,"/ExBV_Schlegel_Inf",&exbv3[0],"a","",exbv3.size());
			ebf::Write(filename1,"/glon",&exbv4[0],"a","degree",exbv4.size());
			ebf::Write(filename1,"/glat",&exbv5[0],"a","degree",exbv5.size());
		}
		else
		{
			cout<<"Skipping Extinction as data already exists"<<endl;
		}
	}

}

void append_radec(const string& filename,const string dir)
{
//	double PI=3.14159;
//	double DTOR=0.0174533;
	ebf::EbfDataInfo dinfo;
	double RADEG=57.2958;
	double PosC[3],PosS[3],alpha,delta;
//	fbvector<float> x(filename+".ebf","r","/Pos3");
	ebf::EbfVector<float> x(filename+".ebf","/px");
	ebf::EbfVector<float> y(filename+".ebf","/py");
	ebf::EbfVector<float> z(filename+".ebf","/pz");
	vector<float> exbv1(x.dim(0),0.0);
	vector<float> exbv2(x.dim(0),0.0);
//	certify(x.dim(1)==3,"Check map dim");
	for(int64_t i=0;i<x.dim(0);++i)
	{
//		cout<<i<<" "<<x.dim(0)<<" "<<x.size()<<" "<<i*x.dim(0)+2<<endl;
//		PosC[0]=x[i*x.dim(1)+0];	PosC[1]=x[i*x.dim(1)+1];	PosC[2]=x[i*x.dim(1)+2];
		PosC[0]=x[i];	PosC[1]=y[i];	PosC[2]=z[i];
		Cot::xyz_to_lbr(PosC,PosS,3);
		PosS[0]*=RADEG;		PosS[1]*=RADEG;
		if(PosS[0]<0)
			PosS[0]=360.0+PosS[0];
		Cot::lb_to_radec(PosS[0],PosS[1],alpha,delta);
		exbv1[i]=alpha;
		exbv2[i]=delta;
	}

	if(x.dim(0) >0)
	{
		std::string filename1=filename+".ebf";
		if(ebf::ContainsKey(filename1,"/ra",dinfo)==0)
		{
			ebf::Write(filename1,"/ra",&exbv1[0],"a","degree",exbv1.size());
			ebf::Write(filename1,"/dec",&exbv2[0],"a","degree",exbv1.size());
		}
		else
		{
			cout<<"Skipping ra,dec addition as data already exists"<<endl;
		}
	}
}


void append1(const string& filename,IsochroneDB& ic,const string& photodir,const string& photosys,const string&magcolorNames)
{
	ebf::EbfDataInfo dinfo;	
	IsoFileDescriptor isofile_info(photodir+"/"+photosys+"/IsoFileDescriptor.txt",photosys,magcolorNames);
	vector<double> x(isofile_info.magnames.size()+5,0.0);
	ebf::EbfVector<float>   age(filename,"/Age");
	ebf::EbfVector<float>   feh(filename,"/FeH");
	ebf::EbfVector<float> smass(filename,"/Smass");

//	for(size_t i=0;i<isofile_info.magnames.size();++i)
//		cout<<i<<" "<<isofile_info.magnames[i]<<endl;

	// initialize fbvectors
	ebf::EbfFile lum,teff,grav,mtip,mact;
	lum.Open(filename+".Atmp0","/lum","w",4,"log solar_luminosity");
	teff.Open(filename+".Atmp1","/teff","w",4,"log kelvin");
	grav.Open(filename+".Atmp2","/grav","w",4,"log 0.01 meter/second/second");
	mact.Open(filename+".Atmp3","/mact","w",4,"solar_mass");
	mtip.Open(filename+".Atmp4","/mtip","w",4,"solar_mass");

	vector<ebf::EbfFile>  mags;
	mags.resize(isofile_info.magnames.size());
	for(size_t k=0;k<isofile_info.magnames.size();++k)
	{
		stringstream sout;
		sout<<filename<<".Atmp"<<k+5;
		mags[k].Open(sout.str(),"/"+photosys+"_"+isofile_info.magnames[k],"w",4,"magnitude");
	}
	// interpolate
	for(size_t i=0;i<age.size();++i)
	{
		ic.interpolateTGM(double(age[i]),double(feh[i]),double(smass[i]),x);
		lum.Write(&x[0],1);
		teff.Write(&x[1],1);
		grav.Write(&x[2],1);
		mact.Write(&x[3],1);
		mtip.Write(&x[4],1);
		for(size_t k=0;k<isofile_info.magnames.size();++k)
		{
			mags[k].Write(&x[k+5],1);
		}
	}
	if(ebf::ContainsKey(filename,"/teff",dinfo)==0)
	{
		lum.SaveTo(filename,"a");
		teff.SaveTo(filename,"a");
		grav.SaveTo(filename,"a");
	}
	else
	{
		cout<<"Skipping"<<" "<<"/teff"<<" "<<"/lum"<<" "<<"/grav"<<endl;
		lum.Remove();
		teff.Remove();
		grav.Remove();
	}
	if(ebf::ContainsKey(filename,"/mact",dinfo)==0)
	{
		mact.SaveTo(filename,"a");
		mtip.SaveTo(filename,"a");
	}
	else
	{
		cout<<"Skipping"<<" "<<"/mact"<<" "<<"/mtip"<<endl;
		mact.Remove();
		mtip.Remove();
	}

	if(ebf::ContainsKey(filename,mags[0].getDataName(),dinfo)==0)
	{
		for(size_t k=0;k<isofile_info.magnames.size();++k)
			mags[k].SaveTo(filename,"a");
	}
	else
	{
		cout<<"Skipping"<<" "<<mags[0].getDataName()<<" ---- "<<mags[isofile_info.magnames.size()-1].getDataName()<<endl;
		for(size_t k=0;k<isofile_info.magnames.size();++k)
			mags[k].Remove();
	}


//	ebf::EbfVector<float>  lum(filename+".Atmp0","/lum","w");
//	ebf::EbfVector<float> teff(filename+".Atmp1","/teff","w");
//	ebf::EbfVector<float> grav(filename+".Atmp2","/grav","w");
//	vector<ebf::EbfVector<float> > mags(isofile_info.magnames.size(),ebf::EbfVector<float>());
//	for(size_t k=0;k<isofile_info.magnames.size();++k)
//	{
//		stringstream sout;
//		sout<<filename<<".Atmp"<<k+3;
//		mags[k].init(sout.str(),"/"+photosys+"_"+isofile_info.magnames[k],"w");
//	}

	// interpolate
//	for(size_t i=0;i<age.size();++i)
//	{
//		ic.interpolateTGM(double(age[i]),double(feh[i]),double(smass[i]),x);
//		lum.push_back(float(x[0]));
//		teff.push_back(float(x[1]));
//		grav.push_back(float(x[2]));
//		for(size_t k=0;k<isofile_info.magnames.size();++k)
//			mags[k].push_back(float(x[k+3]));
//	}
//	teff.flush();	grav.flush();	lum.flush();
//	for(size_t k=0;k<isofile_info.magnames.size();++k)
//		mags[k].flush();

	// copy results to file
//	if(ebf::Ebf_ContainsKey_Info(filename.c_str(),teff.getDataName().c_str(),&dinfo)==0)
//	{
//		ebf::Ebf_Copy(lum.getFileName().c_str(),filename.c_str(),"a","");
//		ebf::Ebf_Copy(teff.getFileName().c_str(),filename.c_str(),"a","");
//		ebf::Ebf_Copy(grav.getFileName().c_str(),filename.c_str(),"a","");
//	}
//	else
//		cout<<"Skipping"<<" "<<teff.getDataName()<<" "<<lum.getDataName()<<" "<<grav.getDataName()<<endl;

//	if(ebf::Ebf_ContainsKey_Info(filename.c_str(),mags[0].getDataName().c_str(),&dinfo)==0)
//	{
//		for(size_t k=0;k<isofile_info.magnames.size();++k)
//			ebf::Ebf_Copy(mags[k].getFileName().c_str(),filename.c_str(),"a","");
//	}
//	else
//		cout<<"Skipping"<<" "<<mags[0].getDataName()<<" ---- "<<mags[isofile_info.magnames.size()-1].getDataName()<<endl;


}




// must be defined outside sutils.cpp
//void getTypeCheck()
//{
//	char x1=0;
//	int x2=0;
//	int64_t x3=0;
//	float x4=0;
//	double x5=0;
//	unsigned int x6=0;
//	int status=0;
//	if((getType(x1)==1)&&(getType(x2)==2)&&(getType(x3)==3)&&(getType(x4)==4)&&(getType(x5)==5)&&(getType(x6)==-1))
//		status=1;
//	if(status==0)
//	{
//		cout<<"function getType not working properly"<<endl;
//		exit(1);
//	}
//}

void runModel(SurveyDesign &sur, IsochroneDB& ic, Parameters &All,
		Interp &vcirc)
{
	if (All.popID != 10)
	{
		Timer timer3;
		int i_min, i_max;
		if (All.popID == -1)
		{
			i_min = 0;
			i_max = 10;
		}
		else
		{
			i_min = All.popID;
			i_max = All.popID + 1;
		}
		//----------------------------------------------------
		for (int i = 0; i < i_min; ++i)
			sur.flush();


		//		cout<<"--------------------------------------------------------"<<endl;
		cout << "Generating populations................" << endl;
		cout << "--------------------------------------------------------"
				<< endl;
		for (int i = i_min; i < i_max; ++i)
//			for (int i = 8; i < 9; ++i)
		{
			timer3.start();
			StellarPopulation sp(i, All.posC, All.warpFlareOn, &vcirc,
					All.option, All.inputDir);
			timer3.print("Time Tree generation/reading =");
			timer3.start();
			sp.spawn(sur, ic, All.fSample);
			cout << left << setw(36) << "Stars spawned = " << setw(12)
					<< sur.nstars << endl;
			sur.flush();
			timer3.print("Time Spawning=");
			cout << "--------------------------------------------------------"
					<< endl;
		}
		for (int i = i_max; i < 10; ++i)
			sur.flush();

	}
	else
	{
//		double ltot=0.0,mtot=0.0;
		Sampler imf(10000, 0.07, 100.0, Chabrier_exp_imf1, 1);
		for (unsigned int i = 0; i < All.sat_list.size(); ++i)
		{
			cout << "------------------------------" << endl;
//			char buf[512];
//			sprintf(buf, "%s%s%02d%s", All.inputDir.c_str(),
//					"StellarHalos/halo", All.sat_list[i].first, "/snapshot.ebf");
//			cout << "Halo No=" << All.sat_list[i].first << "  Sat No="
//					<< All.sat_list[i].second << endl;
//			string fname = buf;
			cout << All.sat_list[i].first << "  Sat No="
					<< All.sat_list[i].second << endl;

			Satellite Sat(All.inputDir+All.sat_list[i].first, 0, All.sat_list[i].second, All.hdim, All.nres);
			Sat.spawn1(sur, imf, ic, All.fSample, All.seed + i + 100);
//			mtot+=Sat.mass();			ltot+=Sat.light();
//			cout << "M/L=" << Sat.mass() / Sat.light() << " " << i << endl;
			cout << "-----------Done---------------" << endl;
		}
		sur.flush();
//		cout<<"Mtot ="<<ltot*2/(All.fSample*mtot)<<" "<<ltot/All.fSample<<" "<<mtot/2.0<<endl;
	}

}

int main(int argc, char **argv)
{
	//	copyright();
	Timer timer1, timer2;
	srand48(13);

//	uint64_t x;
//	x=4101842887655102017;
//	uint64_t y;
//	y=INT64_C(4101842887655102017);
//	cout<<x<<" "<<y<<" "<<INT64_MAX<<endl;
//	exit(1);

	//--------------------
	Parameters All;
	All.setFromArguments(argc, argv);
	//	All.print();
	nrRan = Ran(All.seed+4);
	nrGauss = Normaldev(0.0, 1.0, All.seed);
	Interp vcirc(All.inputDir + "Model/vcirc.dat");

	if (All.option == 0)
	{
		char c;
		cout << "Are you sure you want to create the BHTree file (Y/N): ";
		cin >> c;
		cout << endl;
		if (c == 'Y')
		{
			for (int i = 0; i < 10; ++i)
//			for (int i = 9; i < 10; ++i)
			{
				timer2.start();
				StellarPopulation sp(i, All.posC, All.warpFlareOn, &vcirc,
						All.option, All.inputDir);
				timer2.print("Time Tree generation/reading =");
			}
		}
		else
			cout << "Exiting without making tree" << endl;

	}

	if ((All.fSample > 0.0) && (All.option == 1))
	{
		//----------------------------------------------------
		//		stringstream sout;
		//		sout<<All.outputDir<<All.SuSuffix<<"/"<<All.halosatFile<<".ebf";
		//		string fname=sout.str();
		string fname = All.outputDir + All.SuSuffix + "/" + All.outputFile
				+ ".ebf";
		SurveyDesign sur(fname, All.appMagLimits, All.absMagLimits,
				All.colorLimits, &All);
		sur.setGeometry(All.geometryOption, All.longitude, All.latitude, All.surveyArea);
//		sur.setError(All.ErrorOption, All.sigma_r, All.sigma_vr, All.sigma_mu,All.sigma_fe, All.sigma_al);
		sur.setCenter(All.posC);
		srand48(7);
		//----------------------------------------------------
		timer2.start();
		IsochroneDB ic(All.inputDir + "Isochrones/", "padova/", All.photoSys, All.magcolorNames, 0);
		ic.print();
		timer2.print("Time Isocrhone Reading");
		//-----------------------------------------------------
		timer2.start();


		if(All.fieldTableFile.size()==0)
		{
			runModel(sur,ic,All,vcirc);
		}
		else
		{
			for(size_t i=0;i<All.fieldTable.col[0].size();++i)
			{
				cout<<i<<" long="<< All.fieldTable.col[0][i]<<" lat="<<All.fieldTable.col[1][i]<<endl;
				sur.setGeometry(All.geometryOption, All.fieldTable.col[0][i], All.fieldTable.col[1][i], All.surveyArea,i);
				runModel(sur,ic,All,vcirc);
			}

		}

//		if (All.popID != 10)
//		{
//
//			Timer timer3;
//			int i_min, i_max;
//			if (All.popID == -1)
//			{
//				i_min = 0;
//				i_max = 10;
//			}
//			else
//			{
//				i_min = All.popID;
//				i_max = All.popID + 1;
//			}
//			//----------------------------------------------------
//			for (int i = 0; i < i_min; ++i)
//				sur.flush();
//
//			//		cout<<"--------------------------------------------------------"<<endl;
//			cout << "Generating populations................" << endl;
//			cout << "--------------------------------------------------------"
//					<< endl;
//			for (int i = i_min; i < i_max; ++i)
//			{
//				timer3.start();
//				StellarPopulation sp(i, All.posC, All.warpFlareOn, &vcirc,
//						All.option, All.inputDir);
//				timer3.print("Time Tree generation/reading =");
//				timer3.start();
//				sp.spawn(sur, ic, All.fSample);
//				cout << left << setw(36) << "Stars spawned = " << setw(12)
//						<< sur.nstars << endl;
//				sur.flush();
//				timer3.print("Time Spawning=");
//				cout
//						<< "--------------------------------------------------------"
//						<< endl;
//			}
//			for (int i = i_max; i < 10; ++i)
//				sur.flush();
//			sur.close();
//		}
//		else
//		{
//			Sampler imf(10000,0.07,100.0,Chabrier_exp_imf1,1);
//			for (unsigned int i = 0; i < All.sat_list.size(); ++i)
//			{
//				char buf[512];
//				cout << "------------------------------" << endl;
//				sprintf(buf, "%s%s%02d%s", All.inputDir.c_str(), "StellarHalos/halo",
//						All.sat_list[i].first, "/snapshot.ebf");
//				cout << "Halo No=" << All.sat_list[i].first << "  Sat No="
//						<< All.sat_list[i].second << endl;
//				string fname=buf;
//				Satellite Sat(fname, All.sat_list[i].second,
//						All.sat_list[i].first * 1000 + All.sat_list[i].second,All.hdim,All.nres);
//				Sat.spawn1(sur, imf, ic, All.fSample,All.seed+i+100);
//				cout << "M/L=" << Sat.mass() / Sat.light() << " " << i << endl;
//				cout << "-----------Done---------------" << endl;
//			}
//			sur.flush();
//		}

		sur.close();

	}

	if ((All.option == 2) || (All.option == 1))
	{
		if (All.addstring == "")
		{
			cout << "Calulating magnitudes................" << endl;
			//----------------------------------------------------
			timer2.start();
			IsochroneDB ic(All.inputDir + "Isochrones/", "padova/",
					All.photoSys, All.magcolorNames, 1);
			ic.print();
			timer2.print("Time Isocrhone Reading");
			//----------------------------------------------------
			if (All.option == 1)
				append1(All.outputDir + All.outputFile + ".ebf", ic,
						All.inputDir + "Isochrones/padova/", All.photoSys,
						All.magcolorNames);
			else
				append1(All.outputFile + ".ebf", ic,
						All.inputDir + "Isochrones/padova/", All.photoSys,
						All.magcolorNames);

		}
		else if (All.addstring == "radec")
		{
			append_radec(All.outputFile, All.inputDir);
		}
		//----------------------------------------------------
		//		cout<<"--------------------------------------------------------"<<endl;
	}

	if ((All.option == 3) || (All.option == 1))
	{
		cout << "Calulating Extinction................" << endl;
		timer2.start();
		if (All.option == 1)
			append_ext(All.outputDir + All.outputFile, All.inputDir);
		else
			append_ext(All.outputFile, All.inputDir);
		timer2.print("Time for extinction calculation");
		//		cout<<"--------------------------------------------------------"<<endl;
		//----------------------------------------------------
	}

	timer1.print("Total Time=");

	return 1;
}






