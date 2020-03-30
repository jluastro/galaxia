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

#include "SurveyDesign.h"
#include "Functions.h"

void stringstreamTovchar(stringstream &sout,vector<char> &cv)
{
	int64_t pos = sout.tellg();
	sout.seekg (0, ios::end);
	int length = sout.tellg();
	cv.resize(length+1);
	sout.seekg (0, ios::beg);
	sout.read (&cv[0],length);
	cv[length]=0;
	sout.seekg (pos, ios::beg);
}



SurveyDesign::~SurveyDesign()
{
	if(geo!=NULL)
		delete geo;
// TODO Auto-generated destructor stub
}

void SurveyDesign::setCenter(double *pos)
{
	for (int i=0;i<6;++i)
		posC[i]=pos[i];
}

void SurveyDesign::writeStars(vector<StarParticlef> &Stars1)
{
	for (size_t i = 0; i < Stars1.size(); ++i)
	{
//		ebfA0.push_back(Stars1[i].Pos[0]);
//		ebfA1.push_back(Stars1[i].Pos[1]);
//		ebfA2.push_back(Stars1[i].Pos[2]);
//		ebfA3.push_back(Stars1[i].Pos[3]);
//		ebfA4.push_back(Stars1[i].Pos[4]);
//		ebfA5.push_back(Stars1[i].Pos[5]);
//		ebfA6.push_back(Stars1[i].feh());
//		ebfA7.push_back(Stars1[i].alpha());
//		ebfA8.push_back(Stars1[i].smass());
//		ebfA9.push_back(Stars1[i].age());
//		ebfA10.push_back(Stars1[i].rad());
//		ebfA11.push_back(Stars1[i].mag(0));
//		ebfA12.push_back(Stars1[i].mag(1));
//		ebfA13.push_back(Stars1[i].mag(2));
//		ebfA14.push_back(Stars1[i].popID());
//		ebfA15.push_back(Stars1[i].satID());
//		ebfA16.push_back(fieldID);
//		ebfA17.push_back(Stars1[i].partID());

		ebfA[0].Write(&(Stars1[i].Pos[0]),1);
		ebfA[1].Write(&(Stars1[i].Pos[1]),1);
		ebfA[2].Write(&(Stars1[i].Pos[2]),1);
		ebfA[3].Write(&(Stars1[i].Pos[3]),1);
		ebfA[4].Write(&(Stars1[i].Pos[4]),1);
		ebfA[5].Write(&(Stars1[i].Pos[5]),1);
		ebfA[6].Write(&(Stars1[i].feh()),1);
		ebfA[7].Write(&(Stars1[i].alpha()),1);
		ebfA[8].Write(&(Stars1[i].smass()),1);
		ebfA[9].Write(&(Stars1[i].age()),1);
		ebfA[10].Write(&(Stars1[i].rad()),1);
		ebfA[11].Write(&(Stars1[i].mag(0)),1);
		ebfA[12].Write(&(Stars1[i].mag(1)),1);
		ebfA[13].Write(&(Stars1[i].mag(2)),1);
		ebfA[14].Write(&(Stars1[i].popID()),1);
		ebfA[15].Write(&(Stars1[i].satID()),1);
		ebfA[16].Write(&fieldID,1);
		ebfA[17].Write(&(Stars1[i].partID()),1);

	}
}


int SurveyDesign::push_back(StarParticle &Star)
{

//	ParticleStar pp;
//	pp.Pos = &Star.Pos[0];

	//	if(checkColMag(Star)&&(geo->check1(&pp)))
	if (checkColMag(Star))
	{
//		error->add(Star);
		if (geo->checkP(&Star.pos(0)))
		{
			Star.convert(Starf);
			Stars.push_back(Starf);
			nstars++;
			if(nstars==numeric_limits<int>::max())
			{
				cout<<"ERROR: number of stars larger than maximum int limit- options are"<<endl;
				cout<<"1) Reduce magnitude limits"<<endl;
				cout<<"2) Reduce fSample"<<endl;
				cout<<"3) Make changes in code in file 'SurveyDesign.cpp' to handle larger integers"<<endl;
				exit(1);
			}
			if (Stars.size() == 100000)
			{
//				ebf1.write(&Stars[0], sizeof(Stars[0]) * Stars.size());
				writeStars(Stars);
				Stars.clear();
			}
			return 1;
		}
		else
			return 0;
	}
	return 0;
}

int SurveyDesign::force_push_back(StarParticle &Star)
{
	Star.convert(Starf);
	Stars.push_back(Starf);
	nstars++;
	if(nstars==numeric_limits<int>::max())
	{
		cout<<"ERROR: number of stars larger than maximum int limit- options are"<<endl;
		cout<<"1) Reduce magnitude limits"<<endl;
		cout<<"2) Reduce fSample"<<endl;
		cout<<"3) Make changes in code in file 'SurveyDesign.cpp' to handle larger integers"<<endl;
		exit(1);
	}
	if (Stars.size() == 100000)
	{
		writeStars(Stars);
		Stars.clear();
	}
	return 1;
}

int SurveyDesign::push_check(StarParticle &Star)
{
	if (checkColMag(Star))
	{
		if (geo->checkP(&Star.pos(0)))
			return 1;
		else
			return 0;
	}
	return 0;
}


void SurveyDesign::flush( )
{
//	cout<<Stars.size()<<" "<<nstars<<endl;
	if(Stars.size()>0)
	{
//		ebf1.write(&Stars[0],sizeof(Stars[0])*Stars.size());
		writeStars(Stars);
		Stars.clear();
	}
//	printv(npart.begin(),npart.end());
	npart.push_back(nstars);
//	printv(npart.begin(),npart.end());
	nstars=0;

}

void SurveyDesign::close()
{
	certify(Stars.size()==0);
	nstars=total(npart);
	if(nstars==numeric_limits<int>::max())
	{
		cout<<"ERROR: number of stars larger than maximum int limit- options are"<<endl;
		cout<<"1) Reduce magnitude limits"<<endl;
		cout<<"2) Reduce fSample"<<endl;
		cout<<"3) Make changes in code in file 'SurveyDesign.cpp' to handle larger integers"<<endl;
		exit(1);
	}
//	ebf1.ebfh.dim[0]=nstars;
//	ebf1.close();
//	assert(ebf1.dim(0)==nstars);
//	ebf1.open(outputFile.c_str(),"r+","/Typelist1",2,npart.size());
//	ebf1.sync(&npart[0]);
	cout<<left<<setw(36)<<"Total stars written "<<setw(24)<<nstars<<endl;
	{
		ebf::EbfFile ebf3;

//		for(size_t i=0;i<ebfA.size();++i)
//			ebfA[i].close();

//		ebfcat(ebfA[2].getFilename().c_str(),"/FeH");
//		exit(1);

//		outputFile=outputFileR;		outputFile+=".check";
		stringstream sout;
		sout<<"# File generated by "<<version()<<endl;
//		time_t t1=time(NULL);
//		sout<<"# "<<ctime(&t1);
		sout<<"# <parameterfile>"<<endl;
		vector<char> buf;
		sout<<All->outParameterFile();
		sout<<"# </parameterfile>"<<endl;
		stringstreamTovchar(sout,buf);
		ebf::Write(outputFileR,"/Log",&buf[0],"w","",buf.size());


//		ebf::Write(outputFileR,"/Typelist",&nstars,"a","",1);
//		ebf::Write(outputFileR,"/Typelist1",&npart[0],"a","",npart.size());

//		ebfA0.saveTo(outputFileR,"a");
//		ebfA1.saveTo(outputFileR,"a");
//		ebfA2.saveTo(outputFileR,"a");
//		ebfA3.saveTo(outputFileR,"a");
//		ebfA4.saveTo(outputFileR,"a");
//		ebfA5.saveTo(outputFileR,"a");
//		ebfA6.saveTo(outputFileR,"a");
//		ebfA7.saveTo(outputFileR,"a");
//		ebfA8.saveTo(outputFileR,"a");
//		ebfA9.saveTo(outputFileR,"a");
//		ebfA10.saveTo(outputFileR,"a");
//		ebfA11.saveTo(outputFileR,"a");
//		ebfA12.saveTo(outputFileR,"a");
//		ebfA13.saveTo(outputFileR,"a");
//		ebfA14.saveTo(outputFileR,"a");
//		ebfA15.saveTo(outputFileR,"a");
//		ebfA16.saveTo(outputFileR,"a");
//		ebfA17.saveTo(outputFileR,"a");

		for(size_t i=0;i<ebfA.size();++i)
		{
			ebfA[i].SaveTo(outputFileR,"a");
//			if(printFlag[i]==1)
//				ebfformat2::ebfcopy(ebfA[i].getFilename(),outputFileR,"r+");
//			remove(ebfA[i].getFilename().c_str());
		}

		ebf::Write(outputFileR,"/Center",&posC[0],"a","1000 parsec",6);

		if(All->fieldTable.col.size() != 0)
		{
			ebf::Write(outputFileR,"/Field/Longitude",&(All->fieldTable.col[0][0]),"a","degree",All->fieldTable.col[0].size());
			ebf::Write(outputFileR,"/Field/Latitude",&(All->fieldTable.col[1][0]),"a","degree",All->fieldTable.col[1].size());
		}



//		outputFile=outputFileR;		outputFile+=".tmp";
	}
//	reformat();
//	remove(outputFile.c_str());

	cout<<left<<setw(36)<<"File written- "<<outputFileR<<endl;
//	cout<<"--------------------------------------------------------"<<endl;
}

void SurveyDesign::setError(int option,double sigma_r,double sigma_vr,double sigma_mu,double sigma_fe,double sigma_al)
{
//	error=new SurveyError(option);
//	error->sigma_r=sigma_r;
//	error->sigma_vr=sigma_vr;
//	error->sigma_mu=sigma_mu;
//	error->sigma_fe=sigma_fe;
//	error->sigma_al=sigma_al;
}

void SurveyDesign::setGeometry(int option,double l,double b,double area, int fieldNo)
{
	if(geo!=NULL)
		delete geo;
	fieldID=fieldNo;

	double temp,th,dth;
	double rad2deg=180.0/PI;
	double deg2rad=PI/180.0;
	switch (option)
	{
	case 0:
		geo=new AllSky();
		cout<<left<<setw(36)<<"Using geometry:"<<"All Sky"<<endl;
		break;
	case 1:
		temp=area*(deg2rad*deg2rad)/(2*PI);
		assert((temp>=0.0)&&(temp<=2));
		dth=rad2deg*min(acos(1-temp),PI);
		cout<<left<<setw(36)<<"Using geometry:"<<"Patch at l ,b : ("<<l<<" "<<b<<") d_theta="<<dth<<endl;
		geo=new GPatch(l,b,dth);
		break;
	case 2:
		temp=60.0;
		th=rad2deg*max(acos(cos(temp*deg2rad)+area*(deg2rad*deg2rad)/(2*PI)),0.0);
		dth=temp-th;
		cout<<left<<setw(36)<<"Using geometry:"<<"Angle Cone Strip th dth: "<<th<<" "<<th+dth<<endl;
		geo=new Cone(th,dth);
		break;
	case 3:
		temp=rad2deg*min(acos(area*(3.0/2.0)*(deg2rad*deg2rad)/(2*PI)),PI/2);
		th=temp;
		dth=2*(90.0-th);
		cout<<left<<setw(36)<<"Using geometry:"<<"Angle Wedge X th dth: "<<th<<" "<<th+dth<<endl;
		geo=new Wedge(th,dth,0);
		break;
	case 4:
		temp=rad2deg*min(acos(area*(3.0/2.0)*(deg2rad*deg2rad)/(2*PI)),PI/2);
		th=temp;
		dth=2*(90.0-th);
		cout<<left<<setw(36)<<"Using geometry:"<<"Angle Wedge Y th dth: "<<th<<" "<<th+dth<<endl;
		geo=new Wedge(th,dth,1);
		break;
	case 5:
		th=0.0;
		dth=rad2deg*min(acos(1-area*(deg2rad*deg2rad)/(2*PI)),PI);
		cout<<left<<setw(36)<<"Using geometry:"<<"Angle Cone Patch th dth:"<<th<<" "<<th+dth<<endl;
		geo=new Cone(th,dth);
		break;
	default:
		break;
	}
}

bool SurveyDesign::checkColMag(StarParticle & Star)
{
	if(All->starType==1)
		if((Star.teff()>(4.00279 -0.079*Star.lum()))||(Star.teff()<(4.00279 -0.079*Star.lum()-0.09))||(Star.lum()<1.4)||(Star.lum()>2.0))
			return 0;

	if(All->starType==2)
		if((Star.lum()<1.4)||(Star.lum()>2.0)||(Star.teff()<3.8)||(Star.smass()<0.5)||(Star.smass()>1.0))
			return 0;


	// changed Jan8 2010
//	double dmod=5*log10(Star.rad()*100.0);
//	double dmod=Star.dmod();
//		&&(Star.age()<=age_max)&&(Star.age()>=age_min)

		if((Star.mag(0) >= appMag[0])&&(Star.mag(0) <= appMag[1])&&((Star.mag(1)-Star.mag(2)) >= color[0])&&((Star.mag(1)-Star.mag(2)) <= color[1])&&((Star.mag(0)-Star.dmod()) >= absMag[0])&&((Star.mag(0)-Star.dmod()) <= absMag[1])&&(Star.rad() <=r_max))
		return 1;
	else
		return 0;

//	bool temp=(Star.mag(0) >= appMag[0])&&(Star.mag(0) <= appMag[1]);
//	temp=temp&&((Star.mag(1)-Star.mag(2)) >= color[0])&&((Star.mag(1)-Star.mag(2)) <= color[1]);
//	temp=temp&&((Star.mag(0)-dmod) >= absMag[0])&&((Star.mag(0)-dmod) <= absMag[1]);
//	temp=temp&&(Star.rad() <=r_max);
//	return temp;
}


//bool SurveyDesign::checkAbsMag(double mag1)
//{
//	return (mag1 <= absMag[1])&&(mag1 >= absMag[0]);
//}








