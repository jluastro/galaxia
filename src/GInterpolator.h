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

#ifndef GINTERPOLATOR_H_
#define GINTERPOLATOR_H_
#include "ebfvector.hpp"
#include "Cot.h"

class GInterpolator
{
public:
	GInterpolator():RADEG(57.2958){};
	~GInterpolator();
	void readFromFile(string fname, string tagname)
	{
		string tname = tagname + ".data";
		string tname1;
		ebf::EbfVector<float> dx(fname, tname);
		data.clear();
		data.resize(dx.size(), 0.0);
		for (size_t i = 0; i < dx.size(); ++i)
			data[i] = dx[i];
		dims.clear();
		for (int64_t i = 0; i < dx.rank(); ++i)
			dims.push_back(dx.dim(i));
		key.clear();
		key.resize(dims.size(), 0.0);
		//		ExtBinFormat2 ebf2;
		//		string tname1=tagname+".key0";
		//		ebf2.open(fname.c_str(),"r",tname1.c_str(),0);
		//		if(ebf2.status==1)
		//		{
		//			option=1;
		//			for(size_t i=0;i<dims.size();++i)
		//			{
		//				tname1=tagname+".x"+tostring(i);
		//				fbvector<double> y(fname,"r",tname1);
		//				x[i].resize(y.size());
		//				for(size_t j=0;j<y.size();++j)
		//					x[i][j]=y[j];
		//			}
		//		}
		//		else
		{
			tname1 = tagname + ".xmms";
			ebf::EbfVector<float> y(fname, tname1);
			int k = 0;

			if(dims.size()>1)
			{
			certify(int64_t(dims.size()) == y.dim(0),
					"Dimension mismatch in GInterpolator");
			for (int64_t i = 0; i < y.dim(0); ++i)
			{
				vector<float> v;
				for (int64_t j = 0; j < y.dim(1); ++j)
				{
					v.push_back(y[k]);
					k++;
				}
				xmms.push_back(v);
			}
			}
			else
			{
				certify(int64_t(dims.size()) == 1,
						"Dimension mismatch in GInterpolator");
				for (int64_t i = 0; i < 1; ++i)
				{
					vector<float> v;
					for (int64_t j = 0; j < y.dim(0); ++j)
					{
						v.push_back(y[k]);
						k++;
					}
					xmms.push_back(v);
				}
			}

		}
//		cout<<"Read file: "<<fname<<endl;
	}
	double interpol_xyz(double* Pos)
	{
		double PosS[3];
		Cot::xyz_to_lbr(Pos,PosS,3);
		PosS[0]*=RADEG;		PosS[1]*=RADEG;
		if(PosS[0]<0)
			PosS[0]=360.0+PosS[0];
		PosS[2]=log10(PosS[2]);
//		cout<<PosS[0]<<" "<<PosS[1]<<" "<<PosS[2]<<" "<<endl;
		return interpol(PosS);

	}
	double interpol(const double* Pos)
	{
		int64_t index;
		double value = 0.0;
		for (size_t i = 0; i < dims.size(); ++i)
		{
			key[i] = (Pos[i] - xmms[i][0]) / xmms[i][2];
			if (key[i] < 0)
				key[i] = 0.0;
			if (key[i] >= (dims[i] - 1))
				key[i] = dims[i]-2 + 0.999999;
//			cout<<i<<" "<<key[i]<<" "<<xmms[i][0]<<" "<<xmms[i][1]<<" "<<xmms[i][2]<<" "<<Pos[i]<<endl;
		}

		if (dims.size() == 1)
		{
			for (size_t i = 0; i < 2; ++i)
				{
					index = int(key[0]);
//					certify(index < int64_t(data.size()),"index for data too karge in GInterpolate");
					value = data[index]*(1 - (key[0] -index))+data[index+1]*(key[0] -index);
				}
		}
		else if (dims.size() == 2)
		{
			for (size_t i = 0; i < 2; ++i)
				for (size_t j = 0; j < 2; ++j)
				{
					index = (int(key[0]) + i) * dims[1] + (int(key[1]) + j);
//					certify(index < int64_t(data.size()),"index for data too karge in GInterpolate");
					value += data[index] * fabs(1 - (key[0] - int(key[0]))-i)
							* fabs(1 - (key[1] - int(key[1]))-j);
				}
		}
		else if (dims.size() == 3)
		{
			for (size_t i = 0; i < 2; ++i)
				for (size_t j = 0; j < 2; ++j)
					for (size_t k = 0; k < 2; ++k)
					{
						index = (int(key[0]) + i) * dims[1] * dims[2]
								+ (int(key[1]) + j) * dims[2] + (int(key[2])
								+ k);
//						certify(index < int64_t(data.size()),"index for data too karge in GInterpolate");
						value += data[index] * fabs(1-(key[0] - int(key[0]))-i)
								* fabs(1-(key[1] - int(key[1]))-j) * fabs(1-(key[2] - int(key[2]))-k);
//						cout<<i<<j<<k<<" "<<value<<" "<<data[index]<<" "<<fabs(1-(key[0] - int(key[0]))-i)*fabs(1-(key[1] - int(key[1]))-j)*fabs(1-(key[2] - int(key[2]))-k)<<endl;
					}
		}
		return value;

	}
	void print()
	{
	}
private:
	//	int option;
	double RADEG;
	vector<float> data;
	vector<int64_t> dims;
	vector<vector<float> > xmms;
	vector<double> key;

};

#endif /* GINTERPOLATOR_H_ */
