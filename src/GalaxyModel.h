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



#ifndef GALAXYMODEL_H_
#define GALAXYMODEL_H_

#include "Parameters.h"

class GalaxyModel
{
public:
	GalaxyModel();
	~GalaxyModel();
	void setFromParameterFile(const string fname);
	void checkCompilation( );

	string GalaxiaData;
	double bulge_sigma_r, bulge_sigma_phi, bulge_sigma_z;
	double bulge_x0, bulge_y0, bulge_z0;
	double bulge_alpha, bulge_beta, bulge_gamma;
	double bulge_Rc, bulge_patternspeed;

	std::vector<ParMem>  par_list;
};

#endif /*GALAXYMODEL_H_*/
