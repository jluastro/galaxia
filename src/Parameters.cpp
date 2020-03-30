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

#include "Parameters.h"
#include "CompatibilityChecks.h"

string version()
{
	string s="0.7.1";
	return s;
}


Parameters::Parameters()
{
#define DOUBLE 1
#define STRING 2
#define INT 3


//---General  options--------------------------------------------
    Debug                 =  0;
    SuSuffix="";
    outputFile="";
    nbodyFile="";
    parameterFile="";
    addstring="";

//----General options-----------------------------------------
	par_list.push_back(ParMem("outputFile",&outputFile,STRING,"Examples/galaxy",1));
	par_list.push_back(ParMem("codeDataDir",&inputDir,STRING,CODEDATAPATH,0));
	if(inputDir[inputDir.length()-1]!='/')
		inputDir+="/";
	par_list.push_back(ParMem("outputDir",&outputDir,STRING,"Examples/",1));
	par_list.push_back(ParMem("photoSys",&photoSys,STRING,"UBV",1));
	par_list.push_back(ParMem("magcolorNames",&magcolorNames,STRING,"V,B-V",1));
	par_list.push_back(ParMem("appMagLimits[0]",&appMagLimits[0],DOUBLE,"-100.0",1));
	par_list.push_back(ParMem("appMagLimits[1]",&appMagLimits[1],DOUBLE,"23.0",1));
	par_list.push_back(ParMem("absMagLimits[0]",&absMagLimits[0],DOUBLE,"-100.0",1));
	par_list.push_back(ParMem("absMagLimits[1]",&absMagLimits[1],DOUBLE,"100.0",1));
	par_list.push_back(ParMem("colorLimits[0]",&colorLimits[0],DOUBLE,"-10.0",1));
	par_list.push_back(ParMem("colorLimits[1]",&colorLimits[1],DOUBLE,"10.0",1));
	par_list.push_back(ParMem("geometryOption",&geometryOption,INT,"0",1));
	par_list.push_back(ParMem("starType",&starType,INT,"0",1));    // 1 RR-Lyrae
	par_list.push_back(ParMem("photoError",&photoError,INT,"0",1));  // 2 LSST 3 less error 4 SDSS
	par_list.push_back(ParMem("surveyArea",&surveyArea,DOUBLE,"0.0",1));
	par_list.push_back(ParMem("rSun[0]",&posC[0],DOUBLE,"-8.0",0));
	par_list.push_back(ParMem("rSun[1]",&posC[1],DOUBLE,"0.0",0));
	par_list.push_back(ParMem("rSun[2]",&posC[2],DOUBLE,"0.015",0));
	par_list.push_back(ParMem("vSun[0]",&posC[3],DOUBLE,"11.1",0));
	par_list.push_back(ParMem("vSun[1]",&posC[4],DOUBLE,"239.08",0));  //226.84+12.24
	par_list.push_back(ParMem("vSun[2]",&posC[5],DOUBLE,"7.25",0));

	par_list.push_back(ParMem("fSample",&fSample,DOUBLE,"1.0",1));
	par_list.push_back(ParMem("sigma_r",&sigma_r,DOUBLE,"0.1",0));
	par_list.push_back(ParMem("sigma_vr",&sigma_vr,DOUBLE,"10.0",0));
	par_list.push_back(ParMem("sigma_mu",&sigma_mu,DOUBLE,"200.0",0));
	par_list.push_back(ParMem("sigma_fe",&sigma_fe,DOUBLE,"0.1",0));
	par_list.push_back(ParMem("sigma_al",&sigma_al,DOUBLE,"0.1",0));
	par_list.push_back(ParMem("popID",&popID,INT,"-1",1));
	par_list.push_back(ParMem("warpFlareOn",&warpFlareOn,INT,"1",1));
	par_list.push_back(ParMem("longitude",&longitude,DOUBLE,"0.0",1));
	par_list.push_back(ParMem("latitude",&latitude,DOUBLE,"90.0",1));
	par_list.push_back(ParMem("seed",&seed,INT,"17",1));
	par_list.push_back(ParMem("hdim",&hdim,INT,"6",0));
	par_list.push_back(ParMem("nres",&nres,INT,"64",0));
	par_list.push_back(ParMem("r_max",&r_max,DOUBLE,"1e10",1));


	option=1;
#undef DOUBLE
#undef STRING
#undef INT

}
#undef CODEDATAPATH


Parameters::~Parameters()
{
}

void Parameters::copyright()
{
	cout<<endl<<endl;
	cout<<"--------------------------------------------------------"<<endl;
	cout<<"Galaxia is a code for synthetic modelling of the Milky Way"<<endl;
	cout<<"Coyright (c) 2012 by Sanjib Sharma"<<endl<<endl;
	cout<<"This program is free software: you can redistribute it and/or modify"<<endl;
	cout<<"it under the terms of the GNU Affero General Public License as"<<endl;
	cout<<"published by the Free Software Foundation, either version 3 of the"<<endl;
	cout<<"License, or (at your option) any later version."<<endl<<endl;

	cout<<"This program is distributed in the hope that it will be useful,"<<endl;
	cout<<"but WITHOUT ANY WARRANTY; without even the implied warranty of"<<endl;
	cout<<"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"<<endl;
	cout<<"GNU Affero General Public License for more details."<<endl<<endl;

	cout<<"You should have received a copy of the GNU Affero General Public License"<<endl;
	cout<<"along with this program.  If not, see <http://www.gnu.org/licenses/>."<<endl;
	cout<<"--------------------------------------------------------"<<endl;
	cout<<endl<<endl;
    exit(1);
}

void Parameters::usage(void )
{
	cout<<endl;
	cout<<"NAME:"<<endl;
    cout<<"\t galaxia-"<<version()<<" -a code to generate a synthetic galaxy survey"<<endl;
	cout<<endl;
	cout<<"USAGE:"<<endl;
    cout<<"\t galaxia\t -s parameterfile"<<endl;
    cout<<"\t galaxia\t -r parameterfile"<<endl;
    cout<<"\t galaxia\t -a --psys=photometricSystem filename"<<endl;
    cout<<"\t galaxia\t -r --nfile=haloname [--hdim=3 or 6] parameterfile"<<endl;
    cout<<"\t galaxia\t --copyright "<<endl;
    cout<<"\t galaxia\t --help "<<endl;
	cout<<"DESCRIPTION:"<<endl;
	cout<<"\t -s          "<<"initial setup to generate BHTREE files"<<endl;
	cout<<"\t -r          "<<"run the code to generate stellar data"<<endl;
	cout<<"\t -a          "<<"append catalog file with magnitudes in an alternate photometric system"<<endl;
	cout<<"\t --nfile     "<<"halo02,halo05 etc to sample Bullock Johnston stellar halos"<<endl;
	cout<<"\t --fieldfile "<<"to generate specific fields"<<endl;
	cout<<"\t --hdim      "<<"dimensionality of smoothing lengths, for N-body models only "<<endl;
	cout<<"\t             "<<"6 for with kinematics and 3 for without"<<endl;
	cout<<"\t --copyright "<<"print the copyright and warranty"<<endl;
	 //	cout<<"\t -a --add=radec      "<<"append catalog file with l,b,ra and dec in degrees"<<endl;
	cout<<"CONTACT:"<<endl;
    cout<<"Report bugs to <bugsanjib@gmail.com>."<<endl;
    exit(1);
}



void Parameters::print( )
{
		cout.precision(6);
		cout<<"----------------------------------"<<endl;
		cout<<"Parameter Details"<<endl;
		cout<<"Input Halo Sat File ="<<nbodyFile<<endl;
		cout<<"Photometric System  ="<<photoSys<<endl;
//		cout<<"Error Option        ="<<ErrorOption<<endl;
		cout<<"Geometry Option     ="<<geometryOption<<endl;
		cout<<"Survey Area         ="<<surveyArea<<endl;
		cout<<"Absolute Mag Limits ="<<absMagLimits[0]<<" "<<absMagLimits[1]<<endl;
		cout<<"Apparent Mag Limits ="<<appMagLimits[0]<<" "<<appMagLimits[1]<<endl;
		cout<<"Color Limits        ="<<colorLimits[0]<<" "<<colorLimits[1]<<endl;
		cout<<"fSample             ="<<fSample<<endl;
		//		cout<<left<<setw(12)<<"DensityOn "<<setw(12)<<DensityOn<<endl;
		cout<<"----------------------------------"<<endl;
}


void Parameters::checkCompilation( )
{
	allCompatibilityChecks();

    string s(outputFile);
    if(outputFile.find(".ebf")!=s.npos)
    	outputFile.erase(outputFile.begin()+outputFile.find(".ebf"),outputFile.end());

    if(*(outputDir.rbegin())!='/')
    	outputDir+="/";

    if(*(inputDir.rbegin())!='/')
    	inputDir+="/";

	if(option==0)
	{
		popID=-1;
		fSample=0.0;
	}
	else
	if(option==1)
	{
	}
	else
	if(option==2)
	{
//		if(appendString.size()==0)
//		{
//			cout<<"append string not specified"<<endl;
//			usage();
//		}
//		if(isoFileCustom.size()!=0)
//		{
//			cout<<"Do not specify custom isochrone file"<<endl;
//			usage();
//		}
	}
	else
	if(option==3)
	{
	}
	else
		usage();

	if(nbodyFile.size()!=0)
	{
		outputFile=nbodyFile;
		popID=10;
		load_sat_list();
//		posC[0]=-8.5; posC[1]=0.0; posC[2]=0.0;
//		posC[3]=0.0; posC[4]=0.0; posC[5]=0.0;
	}

	if(fieldTableFile.size()!=0)
		fieldTable.readFromFile(fieldTableFile.c_str());

	if((hdim !=3)and(hdim!=6))
	{
		cout<<"hdim must be 3 or 6"<<endl;
		exit(1);
	}

//	if(sunCentricOn==1)
//	{
//		posC[0]=-8.5; posC[1]=0.0; posC[2]=0.0;
//		posC[3]=0.0; posC[4]=0.0; posC[5]=0.0;
//		cout<<"Sun centering set"<<endl;
//	}
//	else
//	{
//		posC[0]=0.0; posC[1]=0.0; posC[2]=0.0;
//		posC[3]=0.0; posC[4]=0.0; posC[5]=0.0;
//	}
}

//void Parameters::stringToDarr(char *s,T* a,int num)
//{
//	char *s1=NULL;
//	char *s2=NULL;
//	char *c=new char[100];
//	int i=0;
//	s1=s;
//	s2=strchr(s,',')
//	while(s2!=NULL)
//	{
//		strcpy(c,s1,s2-s1);
//		a[i++]=atof(c);
//		s1=s2+1;
//		s2=strchr(s,',')
//	}
//	a[i++]=atof(s1);
//	if(i!=num)
//	{
//		cout<<"Not enough values to fill array"<<endl;
//		exit(1);
//	}
//
//}


void Parameters::load_sat_list()
{
	FILE* fd;
	char buf[512];
	char buf1[512];
	char* cptr;
//	sprintf(buf,"%s%s%s",inputDir,"HaloSatList/",halosatFile);
	string s=inputDir;s+="nbody1/filenames/";s+=nbodyFile;s+=".txt";
	cout<<"Reading Halo Sat File="<<s<<endl;
	sat_list.clear();

	if((fd=fopen(s.c_str(),"r")))
	{
		int sats=-1;
		int x[2];
		string path;
		cptr=fgets(buf,512,fd);
		sscanf(buf,"%s",&buf1[0]);
		path=buf1;
		cptr=fgets(buf,512,fd);
		sscanf(buf,"%d",&sats);
		sscanf(buf,"%d",&x[0]);

//		for(int i=0;i<sats;++i)
		x[1]=0;
		while(1)
		{
			cptr=fgets(buf,512,fd);
			if(feof(fd))
				break;
			if(x[0]==1)
				sscanf(buf,"%s",&buf1[0]);
			else
				sscanf(buf,"%s%d",&buf1[0],&x[1]);
			string s1=path+string(buf1);
			sat_list.push_back(pair<string, int> (s1, x[1]));
//			cout<<"Halo "<<x[0]<<"Sat "<<x[1]<<" "<<feof(fd)<<endl;
		}
		if(int(sat_list.size())!=sats)
		{
			cerr<<"sats!=sat_list.size() check "<<" file"<<endl;
			exit(1);
		}
		fclose(fd);

	}
	if(sat_list.size()==0)
	{
		cerr<<"halo sat list not loaded"<<endl;
		exit(1);
	}

	for(size_t i=0;i<sat_list.size();++i)
	{
		cout<<sat_list[i].first<<" "<<sat_list[i].second<<endl;
	}

	cout<<"No of Satellites   ="<<sat_list.size()<<endl;
}




void Parameters::setFromArguments(int argc, char **argv)
{
//	strcpy(halosatFile, "tnull");
	cout<<"CODEDATAPATH="<<inputDir<<endl;
	option=-1;
	char * c1;
	int i;
		i = 1;
		while (i < argc)
		{
			if ((argv[i][0] == '-') && (argv[i][1] == '-'))
			{

				c1 = strchr(argv[i], '=') + 1;

				if (strcmp(argv[i], "--version") == 0)
				{
					cout << "Version 0.2" << endl;
					exit(1);
				}
				else if (strcmp(argv[i], "--help") == 0)
				{
					usage();
				}
				else if (strcmp(argv[i], "--copyright") == 0)
				{
					copyright();
				}
				else if (strncmp(argv[i], "--nfile=", 8) == 0)
				{
					nbodyFile=c1;
				}
				else if (strncmp(argv[i], "--fieldfile=", 12) == 0)
				{
					fieldTableFile=c1;
				}
				else if (strncmp(argv[i], "--psys=", 7) == 0)
				{
					photoSys=c1;
				}
				else if (strncmp(argv[i], "--hdim=", 6) == 0)
				{
					hdim=atoi(c1);
				}
				else if (strncmp(argv[i], "--nres=", 6) == 0)
				{
					nres=atoi(c1);
				}
				else if (strncmp(argv[i], "--add=", 5) == 0)
				{
					addstring=c1;
				}
//				else if (strncmp(argv[i], "--m_min=", 7) == 0)
//				{
//					m_min=atof(c1);
//				}
				else
					usage();
			}
			else if (argv[i][0] == '-')
			{

				for (size_t j = 1; j < strlen(argv[i]); ++j)
				{
					switch (argv[i][j])
					{
					case 's':
						option=0;
						break;
					case 'r':
						option=1;
						break;
					case 'a':
						option=2;
						break;
					case 'e':
						option=3;
						break;
					default:
						usage();
						break;
					}
				}
			}
			else
//				strcpy(halosatFile, argv[i]);
				parameterFile=argv[i];
			i++;
		}

		if(parameterFile.size()==0)
		{
			cout<<"Parameter file or required argument  not specified"<<endl;
			usage();
		}
		if(outputFile.size()==0)
		{
			cout<<"Output file file not specified"<<endl;
			usage();
		}

		if(option==0)
		{
			if (parameterFile == "nowarp")
				warpFlareOn=0;
			else if (parameterFile == "warp")
				warpFlareOn=1;
			else
			{
				cout<<"required argument should be warp  or nowarp"<<endl;
				usage();
			}

		}
		else if((option==2)||(option==3))
		{
			outputFile=parameterFile;
			parameterFile="";
			magcolorNames="?,?-?";
			if(outputFile.find(".ebf")==string::npos)
			{
				cout<<"The file to be operated on should be a .ebf File"<<endl;
				usage();
			}
		}
		else
			setFromParameterFile(parameterFile);
//		saveParameterFile(halosatFile);



	checkCompilation();

}

void Parameters::saveParameterFile(const string fnameToSave)
{
	ofstream fout(fnameToSave.c_str());
	if(fout.is_open()==1)
	{
		stringstream sout;
		for (size_t i = 0; i < par_list.size(); i++)
		{
			if ((par_list[i].status == 0) && (par_list[i].readOn == 1))
			{
				par_list[i].write(sout);
			}
		}
		fout<<sout.str();
		fout.close();
	}
	else
	{
		cout<<"Cannot open outfile "<<fnameToSave<<endl;
		exit(1);
	}
}


string Parameters::outParameterFile()
{
	stringstream sout;
	for (size_t i = 0; i < par_list.size(); i++)
	{
		if(par_list[i].readOn == 1)
			par_list[i].write(sout);
	}
	for (size_t i = 0; i < par_list.size(); i++)
	{
		if(par_list[i].readOn == 0)
			par_list[i].write(sout);
	}
	return sout.str();
}

void Parameters::setFromParameterFile(const string fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
	FILE *fd;
	char* cptr;
	cout<<left<<setw(36)<<"Reading Parameter file- "<<fname<<endl;
	cout<<"--------------------------------------------------------"<<endl;

	char buf[512], buf1[512], buf2[512], buf3[512];
	int i, j, nt;
	int errorFlag = 0;

	nt = par_list.size();
	if ((fd = fopen(fname.c_str(), "r")))
	{

		{
			while (!feof(fd))
			{
				cptr=fgets(buf, 512, fd);
				if(feof(fd))
					continue;
				if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
					continue;
				if (buf1[0] == '%')
					continue;
				if (buf1[0] == '#')
					continue;
				cout<<left<<setw(24)<<buf1<<" "<<setw(24)<<buf2<<endl;

				for (i = 0, j = -1; i < nt; i++)
					if ((strcmp(buf1, par_list[i].tag.c_str()) == 0)&&(par_list[i].status==0)&&(par_list[i].readOn==1))
					{
						j = i;
						par_list[i].status = 1;
						break;
					}

				if (j >= 0)
				{
					par_list[j].read(buf2);
				}
				else
				{
					cout<<"\nError in file"<<fname<<":   Tag "<<buf1<<" not allowed or multiple defined. "<<j<<endl;
					errorFlag = 1;
				}
			}
		}
		fclose(fd);
//		fclose(fdout);

		//     i= strlen(OutputDir);
		//     if(i>0)
		//       if(OutputDir[i-1]!='/')
		//         strcat(OutputDir, "/");
//		sprintf(buf1, "%s%s", fname, "-usedvalues");
		//    sprintf(buf2, "%s%s", OutputDir, "parameters-usedvalues");
		//    rename(buf1, buf2);
	}
	else
	{
		cout<<"Parameter file "<<fname<<" not found."<<endl;
		errorFlag = 1;
		exit(1);
	}

	for (i = 0; i < nt; i++)
	{
		if ((par_list[i].status==0)&&(par_list[i].readOn==1))
		{
			cout<<"Error. I miss a value for tag '"<<par_list[i].tag<<"' in parameter file '"<<fname<<"'."<<endl;
			errorFlag = 1;
		}
	}

	if (errorFlag)
		exit(1);
	//------------------------------------------------------------

	cout<<"--------------------------------------------------------"<<endl;
//	cout<<"Reading parameterfile Done"<<endl;
#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS

}





//------------------------------------------------------------
