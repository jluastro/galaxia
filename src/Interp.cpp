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

#include "Interp.h"

Interp::Interp()
{
	// TODO Auto-generated constructor stub

}

Interp::~Interp()
{
	// TODO Auto-generated destructor stub
}



void Interp::setFromFile(const string& fname)
{
	   float temp;
	    ifstream fd;
	    int i,no;
	    x.clear();
	    y.clear();
	    fd.open(fname.c_str());
	    if (fd.is_open())
	    {
	    	cout<<left<<setw(36)<<"Reading tabulated values from file- "<<fname<<endl;//" ....."<<flush;
	    	fd>>no;
		for(i=0;i<no;++i)
		{
		    fd>>temp;
		    x.push_back(temp);
		}
		for(i=0;i<no;++i)
		{
		    fd>>temp;
		    y.push_back(temp);
		}
		fd>>temp;

		fd.close();
	    }
	    else
	    {
	    	cout<<"Error opening vcirc file: "<<fname<<endl;
	    	exit(1);
	    }

//	    cout<<"Done"<<endl;
}
