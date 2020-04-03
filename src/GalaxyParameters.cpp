/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Author: Michael Medford                                                   *
 * Copyright (c) 2020 Michael Medford                                        *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of Galaxia. The full Galaxia copyright notice, including*
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and COPYRIGHT, which can be found at the root           *
 * of the source code distribution tree.                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "GalaxyParameters.h"


GalaxyParameters::GalaxyParameters()
{
#define DOUBLE 1
#define STRING 2
#define INT 3
    par_list.push_back(ParMem("GalaxiaData", &GalaxiaData, STRING, "/abc/123", 1));
	par_list.push_back(ParMem("mass", &mass, DOUBLE, "1.2345", 1));
	par_list.push_back(ParMem("ncomp", &ncomp, INT, "0", 1));
#undef DOUBLE
#undef STRING
#undef INT

}
#undef CODEDATAPATH


GalaxyParameters::~GalaxyParameters()
{
}

void GalaxyParameters::print( )
{
		cout.precision(6);
		cout<<"----------------------------------"<<endl;
		cout<<"GalaxyParameters"<<endl;
		cout<<"GalaxiaData = "<<GalaxiaData<<endl;
		cout<<"mass = "<<mass<<endl;
		cout<<"ncomp = "<<ncomp<<endl;
		cout<<"----------------------------------"<<endl;
}


void GalaxyParameters::checkCompilation( )
{
    if(*(GalaxiaData.rbegin())!='/')
    	GalaxiaData+="/";
}

void GalaxyParameters::setFromParameterFile(const string fname)
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
					cout<<"\nError in file "<<fname<<":   Tag "<<buf1<<" not allowed or multiple defined. "<<j<<endl;
					errorFlag = 1;
				}
			}
		}
		fclose(fd);
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

	cout<<"--------------------------------------------------------"<<endl;
	checkCompilation();

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
}
