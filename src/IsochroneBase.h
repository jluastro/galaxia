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

#ifndef ISOCHRONEBASE_H_
#define ISOCHRONEBASE_H_
#include "Isochrone.h"

class IsoFileDescriptor
{
public:
	IsoFileDescriptor(const string &fname,const string& photoSys,const string &magcolorNames);
	void setExtraFields(int extraFieldsOn)
	{
		if((extraFieldsOn==0)||(extraFieldsOn==1))
		{
			if(extraFieldsOn==1)
			for(size_t i=0;i<magnames.size();++i)
				extraid.push_back(7+i);
		}
		else
		{
			cout<<"extraFieldsOn must be 0 1"<<endl;
			exit(1);
		}
	}
	int fields;
	int magid[3];
	vector<int> extraid;
	vector<string> magnames;
};

class icNeighbor
{
public:
//	double x[4];
	double prob;
	Isochrone* ic;
};


class IsochroneBase
{
public:
	IsochroneBase(const string& inputDir,const string& dirname,const string& photoSys,int extraFieldsOn,const string& magcolorNames);
	virtual ~IsochroneBase();
	void print();
	int magswap;
private:
	int readfile(const string& fname,double alpha1,double feH1,const IsoFileDescriptor &iso_fileinfo,int dwarfOn);
protected:
	void readIsocrhones(const string& inputDir,const string& dirname,const string& photoSys,int extraFieldsOn,const string& magcolorNames);
	inline int index(int i,int j,int k)
	{
	//	return ((k*FeH.size()+j)*Age.size()+i);
		return (j*Age.size()+i);
	}
	inline void reverse_index(int &i,int &j,int &k,int ic_index)
	{
		i=ic_index%int(Age.size());
		j=ic_index/int(Age.size());
		k=0;
	}
	vector<double> Age,FeH,Alpha;
//	double dAge;//,dFeH,dAlpha;
	vector<int>	im;
	vector<Isochrone>	icv;
	vector<icNeighbor> icQ;
};

#endif /* ISOCHRONEBASE_H_ */
