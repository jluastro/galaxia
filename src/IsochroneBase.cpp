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

#include<algorithm>
#include "IsochroneBase.h"
#include "TableInterpolator.h"
#include "ebf.hpp"



//IsoFileDescriptor::IsoFileDescriptor(const string fname,const string& magcolorNames)
//{
//
//	ifstream fd;
//	fd.open(fname.c_str());
//	if (fd.is_open())
//	{
//		fd>>fields>>magid[0]>>magid[1]>>magid[2];
//		int nmags1;
//		fd>>nmags1;
//		string s;
//		for(int i=0;i<nmags1;++i)
//		{
//			fd>>s;
//			magnames.push_back(s);
////			cout<<i<<" "<<s<<endl;
//		}
//
//	    if(fd.good()!=1)
//		{
//			cout << "fscanf problem in IsoFileDescriptor: Incorrect no of params" << endl;
//			exit(1);
//		}
//	    fd.close();
//
//	    vector<string> sv;
//	    stringSplit(magcolorNames," ,-",sv);
//	    certify(sv.size()==3,"mag color string not correct");
//	    cout<<sv[0]<<" "<<sv[1]<<" "<<sv[2]<<endl;
//	    for(size_t i=0;i<sv.size();++i)
//	    {
//	    	vector<string>::iterator it=find(magnames.begin(),magnames.end(),sv[i]);
//	    	certify(it!=magnames.end(),"color or magnitude name not in list");
//	    	magid[i]=7+it-magnames.begin();
//	    }
//    	cout<<sv[0]<<"("<<magid[0]<<") "<<sv[1]<<"("<<magid[1]<<") "<<sv[2]<<"("<<magid[2]<<")"<<endl;
//
//		int fields_r = *max_element(magid, magid + 3);
//		if ((fields_r < 0) || (fields_r >= fields) || (fields_r >= 16))
//		{
//			cout << "Something wrong with field descriptor " << fields << " "
//					<< fields_r << endl;
//			exit(1);
//		}
//	}
//	else
//	{
//		cout << "Isocrhone File Descriptor \"" << fname << "\" not found" << endl;
//		exit(1);
//	}
//
//
//}


IsoFileDescriptor::IsoFileDescriptor(const string &fname,const string& photoSys,const string& magcolorNames)
{

	ifstream fd;
//	cout<<fname<<endl;
	fd.open(fname.c_str());
	if (fd.is_open())
	{
		string photoSys1,s;
		int startid,nmags1;

		while(fd.good()==1)
		{
			getline(fd,s);
			stringstream ss(s);
			ss>>photoSys1>>fields>>startid>>nmags1;
			{
				vector<string> sv1;
				stringSplit(s," ",sv1);
				if(int(sv1.size()) != (nmags1+4))
				{
					cout << "scan problem in IsoFileDescriptor: Incorrect no of params" << endl;
					exit(1);
				}
			}
			magnames.clear();
			for(int i=0;i<nmags1;++i)
			{
				ss>>s;
				magnames.push_back(s);
//							cout<<i<<" "<<s<<" "<<ss.good()<<" "<<ss.eof()<<endl;
			}
//			cout<<photoSys1<<" "<<fields<<" "<<startid<<" "<<nmags1<<" "<<magnames.size()<<" "<<ss.good()<<endl;
			if(photoSys1==photoSys)
				break;
		}

	    if(fd.good()!=1)
		{
			cout << "scan problem in IsoFileDescriptor: Incorrect no of params" << endl;
			exit(1);
		}
	    fd.close();

	    vector<string> sv;
	    stringSplit(magcolorNames," ,-",sv);
	    certify(sv.size()==3,"mag color string not correct");
//	    cout<<sv[0]<<" "<<sv[1]<<" "<<sv[2]<<endl;

	    if(sv[0]=="?")
	    	sv[0]=magnames[0];
	    if(sv[1]=="?")
	    	sv[1]=magnames[1];
	    if(sv[2]=="?")
	    	sv[2]=magnames[0];


	    for(size_t i=0;i<sv.size();++i)
	    {
	    	vector<string>::iterator it=find(magnames.begin(),magnames.end(),sv[i]);
	    	certify(it!=magnames.end(),"color or magnitude name not in list");
	    	magid[i]=7+it-magnames.begin();
	    }
//    	cout<<sv[0]<<"("<<magid[0]<<"), "<<sv[1]<<"("<<magid[1]<<")-"<<sv[2]<<"("<<magid[2]<<")"<<endl;

		int fields_r = *max_element(magid, magid + 3);
		if ((fields_r < 0) || (fields_r >= fields) || (fields_r >= 16))
		{
			cout << "Something wrong with field descriptor " << fields << " "
					<< fields_r << endl;
			exit(1);
		}
	}
	else
	{
		cout << "Isocrhone File Descriptor \"" << fname << "\" not found" << endl;
		exit(1);
	}


}



void IsochroneBase::print()
{
	cout<<left<<setw(36)<<"Isochrone Grid Size:"<<"(Age bins="<<Age.size()<<",Feh bins="<<FeH.size()<<",Alpha bins="<<Alpha.size()<<")"<<endl;
//	printv(Age,"Age");
//	printv(FeH,"FeH");
//	printv(Alpha,"Alpha");
}




IsochroneBase::IsochroneBase(const string& inputDir,const string& dirname,const string& photoSys,int extraFieldsOn,const string &magcolorNames)
{
	magswap=0;
	im.resize(300);
	for(int i=0;i<300;++i)
		im[i]=i;
	readIsocrhones(inputDir,dirname,photoSys,extraFieldsOn,magcolorNames);
}


IsochroneBase::~IsochroneBase()
{
	// TODO Auto-generated destructor stub
}

void IsochroneBase:: readIsocrhones(const string& inputDir,const string& dirname,const string& photoSys,int extraFieldsOn,const string &magcolorNames)
{
	int dwarfOn=1;
	cout<<left<<setw(36)<<"Reading Isochrones from dir- "<<inputDir+dirname<<photoSys<<endl;
	if (ebf::ebfutils::FileExists(inputDir+"BolCorr/"+photoSys+"/interp_keys.txt") == 0)
	{
		dwarfOn=0;
	}

	double  al=0.0;
//	double fe[]={0.0001,0.001,0.01,0.02,0.03};
//	double fe1[]={0.0001,0.0002,0.0005,0.0007,0.0009,0.0012,0.0016,0.002,0.0024,0.003,0.004,0.006,0.008,0.01,0.014,0.018,0.024,0.03};

	double fe1[]={0.0001,0.0002,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.0012,0.0014,0.0016,0.0018,0.002,0.0022,0.0024,0.0026,0.003,0.0034,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.012,0.014,0.016,0.018,0.02,0.024,0.028,0.03};
	vector<double> fe(fe1,fe1+sizeof(fe1)/sizeof(fe1[0]));
	for(size_t i=0;i<fe.size();++i)
	{
		FeH.push_back(log10(fe[i]/0.019));
	}

//	vector<double> fe(64);
//	for(size_t i=0;i<fe.size();++i)
//	{
//		fe[i]=pow(10.0,-4+0.0393194*i*64/fe.size());
//		stringstream ss;
//		ss<<setprecision(5)<<setiosflags(ios_base::fixed)<<fe[i];
//		ss>>fe[i];
//		FeH.push_back(log10(fe[i]/0.019));
////		FeH.push_back(log10(1e-4/0.019)+i*0.0393194);
//	}

	Alpha.push_back(0.0);
	string buf=inputDir+dirname+"/"+photoSys+"/IsoFileDescriptor.txt";
//	string buf;
//	buf=inputDir+dirname;buf+="/IsoFileDescriptor.txt";


	IsoFileDescriptor iso_fileinfo(buf,photoSys,magcolorNames);
	iso_fileinfo.setExtraFields(extraFieldsOn);

	if(iso_fileinfo.magid[0]==iso_fileinfo.magid[1])
		magswap=1;
	if(iso_fileinfo.magid[0]==iso_fileinfo.magid[2])
		magswap=2;

	if(extraFieldsOn>0)
	{
		certify(iso_fileinfo.extraid.size()==(iso_fileinfo.magnames.size()),"extraid vs nmags match");
	}

	int n_ages=0;
	for(unsigned int i=0;i<FeH.size();i++)
	{
		stringstream sout;
		sout<<inputDir<<dirname<<photoSys<<"/output_"<<fixed<<setprecision(6)<<fe[i]<<".dat";
		if(i==0)
			n_ages=readfile(sout.str(),al,fe[i],iso_fileinfo,dwarfOn);
		else
		{
			assert(n_ages==readfile(sout.str(),al,fe[i],iso_fileinfo,dwarfOn));
		}
	}


	if (dwarfOn==1)
	{
	//---------------------------------------------------------
	//  compute magnitude of dwarfs
	string dir1=inputDir+"BolCorr/"+photoSys+"/";
	TableInterpolatorBC  tbc(dir1);
	string dir2=inputDir+"Chabrier/";
	TableInterpolator  tb(dir2);
	double Pos1[6];
	for(size_t i=0;i<icv.size();++i)
	{
		for(size_t j=0;j<icv[i].m.size();++j)
		{
			if(icv[i].m[j]<0.15)
			{
//				cout<<j<<" b "<<icv[i].age<<" "<<icv[i].FeH<<" "<<icv[i].m[j]<<" "<<icv[i].Teff[j]<<" "<<icv[i].B[j]<<" "<<icv[i].V[j]<<" "<<icv[i].I[j]<<endl;
				Pos1[0]=icv[i].age;
				Pos1[1]=icv[i].m[j];
				tb.interpol(Pos1);
				icv[i].Teff[j]=log10(tb.value[0]);
				icv[i].Lum[j] =tb.value[1];
				icv[i].Grav[j]=tb.value[2];

				Pos1[0]=icv[i].FeH;
				Pos1[1]=pow(10.0,icv[i].Teff[j]);
				Pos1[2]=icv[i].Grav[j];
				tbc.interpol(Pos1);
				icv[i].Mag0[j]=4.77-2.5*icv[i].Lum[j]-tbc.value[iso_fileinfo.magid[0]-7];
				icv[i].Mag1[j]=4.77-2.5*icv[i].Lum[j]-tbc.value[iso_fileinfo.magid[1]-7];
				icv[i].Mag2[j]=4.77-2.5*icv[i].Lum[j]-tbc.value[iso_fileinfo.magid[2]-7];


				for(size_t k=0;k<iso_fileinfo.extraid.size();++k)
					icv[i].Mags[k][j]=4.77-2.5*icv[i].Lum[j]-tbc.value[iso_fileinfo.extraid[k]-7];

			}
			else
				break;
		}
	}
	}
	//---------------------------------------------------------
	// changed from with readfile (last one might be missing)
	for(size_t i=0;i<icv.size();++i)
	{
		icv[i].setMon();
		icv[i].setTip();
//		cout<<icv[i].m[0]<<" "<<icv[i].Vmon[0]<<endl;
	}

// few checks
	assert(icv.size()==n_ages*FeH.size());
	for(size_t i=1;i<FeH.size();i++)
	{
		assert(icv[0].age==icv[i*n_ages].age);
		assert(icv[n_ages-1].age==icv[(i+1)*n_ages-1].age);
	}
// set Age vector
	assert(icv[0].FeH==icv[1].FeH);
//	assert(fabs(icv[0].age-6.6)<1e-5);  // just making sure not necessary
	double temp=icv[1].age-icv[0].age;
	for(int i=0;i<n_ages;i++)
		Age.push_back(icv[0].age+i*temp);


//	Age.push_back(6.6+i*temp);
//	dAge=0.02;
//	dAlpha=0.0;
}

int IsochroneBase:: readfile(const string& fname,double alpha1,double feH1,const IsoFileDescriptor &iso_fileinfo,int dwarfOn)
{
	char buf[512];
	char* cptr;
	FILE* fd=NULL;
	char buf_check[512];
	vector<float> x(iso_fileinfo.fields,0.0);
	certify(int(x.size())>=16,"number of fields greater than equal to 16");
	certify(int(iso_fileinfo.extraid.size())<=9,"number of magnames grater than 16");
	Isochrone icData;
	icData.age=-1.0;
	icData.FeH=log10(feH1/0.019);
	icData.alpha=alpha1;
//	icData.addDwarfs(int(iso_fileinfo.extraid.size()));
	icData.addDwarfs(int(iso_fileinfo.extraid.size()),dwarfOn);


	int i=0;
	if((fd=fopen(fname.c_str(),"r")))
	{
//		while(!feof(fd))
		while(1)
		{
			cptr=fgets(buf,512,fd);
			if(feof(fd))
				break;
//			cout<<"check "<<strlen(buf)<<endl;
			if(buf[0]=='#')
			    continue;
//			if(iso_fileinfo.fields==21)
//			sscanf(buf,"%G%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&x[6],&x[7],&x[8],&x[9],&x[10],&x[11],&x[12],&x[13],&x[14],&x[15],&x[16],&x[17],&x[18],&x[19],&x[20]);
//			if(iso_fileinfo.fields==16)
			sscanf(buf,"%G%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5],&x[6],&x[7],&x[8],&x[9],&x[10],&x[11],&x[12],&x[13],&x[14],&x[15]);

			sprintf(buf_check,"%f ",x[2]);
			if(strncmp(buf_check,"nan",3)==0)
				continue;
			if(double(x[0])!=icData.age)
			{
				icData.age=x[0];
				icv.push_back(icData);
				i++;
			}

			icv.back().m.push_back(x[1]);
			icv.back().Mact.push_back(x[2]);
			icv.back().Mag0.push_back(x[iso_fileinfo.magid[0]]);
			icv.back().Mag1.push_back(x[iso_fileinfo.magid[1]]);
			icv.back().Mag2.push_back(x[iso_fileinfo.magid[2]]);
			icv.back().Lum.push_back(x[3]);
			icv.back().Teff.push_back(x[4]);
			icv.back().Grav.push_back(x[5]);

			for(size_t k=0;k<iso_fileinfo.extraid.size();++k)
			{
				icv.back().Mags[k].push_back(x[iso_fileinfo.extraid[k]]);
			}

		}
		fclose(fd);
	}
	else
	{
		cout<<"file not found "<<fname<<endl;
		exit(1);
	}
	return i;
}


