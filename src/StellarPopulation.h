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

#ifndef STELLARPOPULATION_H_
#define STELLARPOPULATION_H_
#include"BHTree.h"
#include"IsochroneDB.h"
#include"SurveyDesign.h"

class StellarPopulation
{
public:
	StellarPopulation(int i,const double* posC,int warpFlareOn1,Interp *vcircP1,int option,const string &inputDir);
	void spawn(SurveyDesign &sur,IsochroneDB &ic,double fSample);
	virtual ~StellarPopulation();
private:
	StarParticle Star;
	BHTree BHT;
	Population* cpop;
};

#endif /* STELLARPOPULATION_H_ */
