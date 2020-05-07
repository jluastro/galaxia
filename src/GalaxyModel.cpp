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

#include "GalaxyModel.h"


GalaxyModel::GalaxyModel()
{
#define DOUBLE 1
#define STRING 2
#define INT 3
    par_list.push_back(ParMem("GalaxiaData", &GalaxiaData, STRING, "", 1));
	par_list.push_back(ParMem("bulge_sigma_r", &bulge_sigma_r, DOUBLE, "0", 1));
	par_list.push_back(ParMem("bulge_sigma_phi", &bulge_sigma_phi, DOUBLE, "0", 1));
	par_list.push_back(ParMem("bulge_sigma_z", &bulge_sigma_z, DOUBLE, "0", 1));
	par_list.push_back(ParMem("bulge_x0", &bulge_x0, DOUBLE, "0", 1));
	par_list.push_back(ParMem("bulge_y0", &bulge_y0, DOUBLE, "0", 1));
	par_list.push_back(ParMem("bulge_z0", &bulge_z0, DOUBLE, "0", 1));
	par_list.push_back(ParMem("bulge_alpha", &bulge_alpha, DOUBLE, "0", 1));
	par_list.push_back(ParMem("bulge_beta", &bulge_beta, DOUBLE, "0", 1));
	par_list.push_back(ParMem("bulge_gamma", &bulge_gamma, DOUBLE, "0", 1));
	par_list.push_back(ParMem("bulge_Rc", &bulge_Rc, DOUBLE, "0", 1));
	par_list.push_back(ParMem("bulge_patternspeed", &bulge_patternspeed, DOUBLE, "0", 1));
#undef DOUBLE
#undef STRING
#undef INT
}



GalaxyModel::~GalaxyModel()
{
}

void GalaxyModel::checkCompilation( )
{
    if(*(GalaxiaData.rbegin())!='/')
    	GalaxiaData+="/";
}

void GalaxyModel::setFromParameterFile(const string fname)
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
