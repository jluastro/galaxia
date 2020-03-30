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

#ifndef SURVEYDESIGN_H_
#define SURVEYDESIGN_H_
#include "Geometry.h"
#include "Parameters.h"
#include "SurveyError.h"
#include "ebfvector.hpp"


class SurveyDesign
{
public:

	SurveyDesign(const string fname, double appMagLimits[2],double absMagLimits[2],double colorLimits[2],Parameters* All1)
	{
		All=All1;
//		sunCentricOn=sunCentricOn1;
//		colorMagCutOption=colorMagCutOption1;
		appMag[0]=appMagLimits[0];
		appMag[1]=appMagLimits[1];
		absMag[0]=absMagLimits[0];
		absMag[1]=absMagLimits[1];
		color[0]=colorLimits[0];
		color[1]=colorLimits[1];
		r_max=All->r_max;
//		age_min=All->age_min;
//		age_max=All->age_max;
//		if(sunCentricOn==1)
//		{
//			posC[0]=-8.5; posC[1]=0.0; posC[2]=0.0;
//			posC[3]=0.0; posC[4]=0.0; posC[5]=0.0;
//			cout<<"Sun centering set"<<endl;
//		}
//		else
//		{
//			posC[0]=0.0; posC[1]=0.0; posC[2]=0.0;
//			posC[3]=0.0; posC[4]=0.0; posC[5]=0.0;
//		}
		outputFileR=fname;
		outputFile=outputFileR;		outputFile+=".tmp";
//		sprintf(outputFile,"%s%s",outputFileR,".tmp");
//		strcpy(outputFile,fname);
//		sprintf(outputFile,"%s%s",outputFileR,".tmp");
//		StarParticlef Starf;
//		ebf1.open(outputFile.c_str(),"w","/Pos",4,0,sizeof(Starf)/sizeof(Starf.Pos[0]));

		printFlag.clear();
		printFlag.resize(18,1);

//		 ebfA0.init(outputFileR+".tmp0","/px","w");
//		 ebfA1.init(outputFileR+".tmp1","/py","w");
//		 ebfA2.init(outputFileR+".tmp2","/pz","w");
//		 ebfA3.init(outputFileR+".tmp3","/vx","w");
//		 ebfA4.init(outputFileR+".tmp4","/vy","w");
//		 ebfA5.init(outputFileR+".tmp5","/vz","w");
//		 ebfA6.init(outputFileR+".tmp6","/feh","w");
//		 ebfA7.init(outputFileR+".tmp7","/alpha","w");
//		 ebfA8.init(outputFileR+".tmp8","/smass","w");
//		 ebfA9.init(outputFileR+".tmp9","/age","w");
//		 ebfA10.init(outputFileR+".tmp10","/rad","w");
//		 ebfA11.init(outputFileR+".tmp11","/mag0","w");
//		 ebfA12.init(outputFileR+".tmp12","/mag1","w");
//		 ebfA13.init(outputFileR+".tmp13","/mag2","w");
//		 ebfA14.init(outputFileR+".tmp14","/popID","w");
//		 ebfA15.init(outputFileR+".tmp15","/satID","w");
//		 ebfA16.init(outputFileR+".tmp16","/fieldID","w");
//		 ebfA17.init(outputFileR+".tmp17","/partID","w");


		ebfA.clear();
		ebfA.resize(18);
		ebfA[0].Open(outputFileR+".tmp0","/px","w",4,"1000 parsec");
		ebfA[1].Open(outputFileR+".tmp1","/py","w",4,"1000 parsec");
		ebfA[2].Open(outputFileR+".tmp2","/pz","w",4,"1000 parsec");
		ebfA[3].Open(outputFileR+".tmp3","/vx","w",4,"1000 meter/second");
		ebfA[4].Open(outputFileR+".tmp4","/vy","w",4,"1000 meter/second");
		ebfA[5].Open(outputFileR+".tmp5","/vz","w",4,"1000 meter/second");
		ebfA[6].Open(outputFileR+".tmp6","/feh","w",4);
		ebfA[7].Open(outputFileR+".tmp7","/alpha","w",4);
		ebfA[8].Open(outputFileR+".tmp8","/smass","w",4,"solar_mass");
		ebfA[9].Open(outputFileR+".tmp9","/age","w",4,"log year");
		ebfA[10].Open(outputFileR+".tmp10","/rad","w",4,"1000 parsec");
		ebfA[11].Open(outputFileR+".tmp11","/mag0","w",4);
		ebfA[12].Open(outputFileR+".tmp12","/mag1","w",4);
		ebfA[13].Open(outputFileR+".tmp13","/mag2","w",4);
		ebfA[14].Open(outputFileR+".tmp14","/popid","w",2);
		ebfA[15].Open(outputFileR+".tmp15","/satid","w",2);
		ebfA[16].Open(outputFileR+".tmp16","/fieldid","w",2);
		ebfA[17].Open(outputFileR+".tmp17","/partid","w",2);


		if(All->fieldTable.col.size() == 0)
			printFlag[12]=0;

		outputFile=outputFileR;		outputFile+=".tmp";

		Stars.reserve(100000);
		nstars=0;
		fieldID=0;
		geo=NULL;
	}
	~SurveyDesign();
	Geometry* geo;
//	SurveyError* error;
//	Sampler imf;
//	Interpolator ic;
	int force_push_back(StarParticle &Star);
	int push_back(StarParticle &Star);
	int push_check(StarParticle &Star);
	void flush();
	void close();
	bool checkColMag(StarParticle &Star);
//	bool checkAbsMag(double mag1);
//	void reformat();
//	void reformat1();
//	void reformat2();
	void writeStars(vector<StarParticlef> &Stars1);
	string outputFileR;
	string outputFile;
//	ExtBinFormat ebf1;
	vector<ebf::EbfFile> ebfA;
//	ebf::EbfVector<float> ebfA0,ebfA1,ebfA2,ebfA3,ebfA4,ebfA5,ebfA6,ebfA7,ebfA8,ebfA9,ebfA10,ebfA11,ebfA12,ebfA13;
//	ebf::EbfVector<int> ebfA14,ebfA15,ebfA16,ebfA17;
	vector<int> printFlag;
	double appMag[2];
	double absMag[2];
	double color[2];
	double r_max;//,age_min,age_max;
	Parameters* All;
//	int sunCentricOn,colorMagCutOption;
	double posC[6];
//	double velC[3];
	vector<StarParticlef> Stars;
	int nstars;
	StarParticlef Starf;
	vector<int> npart,npartc;
	int fieldID;
	void setGeometry(int option,double l,double b,double area,int fieldNo=0);
	void setError(int option,double sigma_r,double sigma_vr,double sigma_mu,double sigma_fe,double sigma_al);
	void setCenter(double *pos);
};



#endif /* SURVEYDESIGN_H_ */
