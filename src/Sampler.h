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

#ifndef SAMPLER_H_
#define SAMPLER_H_

#include "sutils.h"

class Sampler
{
public:
//	Sampler(int nsize,double xmin,double xmax,double (* func) (double x));
	Sampler(int nsize,double xmin,double xmax,double (* func) (double x),int optionlinlog);
	Sampler(vector<double> &x,vector<double> &y);
	Sampler(vector<double> &x);
	Sampler(const string fname);
	~Sampler();
	void normalize();
	void calculateCpd( );
	double rand( );
	void print( );
	void plot( );
	vector<double> randv(int nsize1);
	void setRange(double xmin1,double xmax1);
	void setSeed(int64_t seed);
	void getFacv(vector<double> &x1_a,vector<double> &x2_a,vector<double> &fac_a);
	double getFac(double xmin1,double xmax1);
	//	double p[3];
//	void linspace(vector<double>&x,double xmin,double xmax);
//	void logspace(vector<double> &x,double xmin,double xmax);
//	double int_tabulated(vector<double> &x,vector<double> &y,double xmin,double xmax);
//	double int_tabulated(vector<double> &x,vector<double> &y);
	double (*function) (double x);
//	double par_per_lum;
	double meanx,x_min,x_max;
	double cpd_min,cpd_max;
	std::vector<double> xd;
	std::vector<double> x;
	std::vector<double> px;
	std::vector<double> cpd;
};

#endif /*SAMPLER_H_*/
