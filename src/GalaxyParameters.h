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



#ifndef GALAXYPARAMETERS_H_
#define GALAXYPARAMETERS_H_

#include "Parameters.h"

class GalaxyParameters
{
public:
	GalaxyParameters();
	~GalaxyParameters();
	void setFromParameterFile(const string fname);
	void print( );
	void checkCompilation( );
	string GalaxiaData;
	double mass;
	int ncomp;
	std::vector<ParMem>  par_list;
};

#endif /*GALAXYPARAMETERS_H_*/
