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

#ifndef TABLEINTERPOLATOR_H_
#define TABLEINTERPOLATOR_H_

#include<algorithm>
#include "Interp.h"
#include "Functions.h"

class Neighbor
{
public:
	double prob;
	int ig;
};


class GridTableInterpolatorBase
{
public:
	GridTableInterpolatorBase(string fname);
	virtual ~GridTableInterpolatorBase();
	virtual void interpol(const double* Pos)=0;
	void printBase();
	inline int index(int* ijk)
	{
		int temp=ijk[keys.size()-1];
		for(int k=int(keys.size()-2);k>=0;--k)
			temp=(temp*keys[k-1].size()+k);
		return temp;
	}
	inline void reverse_index(int* ijk,int ic_index)
	{
		ijk[0]=ic_index;
		for(int k=1;k<int(keys.size())-1;++k)
		{
			ijk[k]=ijk[k-1]/int(keys[k-1].size());
			ijk[k-1]=ijk[k-1]%int(keys[k-1].size());
		}
	}
	void setNeighbors(const double* Pos);
	vector<double> value;
protected:
	vector<int> iv;
	vector<vector<double> > keys;
	vector<double> dkeys,cPos,du;
	vector<string> file_list;
	vector<Neighbor> nv;
};

class AsciiTable
{
public:
	AsciiTable(const char* fname);
	bool check_comment();
	void next();
	~AsciiTable();
	string _name;
	FILE* fd;
	char buf[1024];
	int rows,columns,dim1,ic,scan_status,massid;
	vector<int> dimid;
	vector<double> x;
	vector<float> tempv;
};


class Table
{
public:
	Table(){};
	void initialize(int rows, int columns)
	{
		if(col.size()>0)
		{
			for(size_t k=0;k<col.size();++k)
				col[k].clear();
			col.clear();
		}
		for(int k=0;k<columns;++k)
			col.push_back(vector<double>());
		for(int k=0;k<columns;++k)
			col[k].resize(rows,0.0);

	}
	void readFromFile(const char* fname)
	{
		AsciiTable at(fname);
		initialize(at.rows,at.columns);
		for(int i=0;i<at.rows;++i)
		{
			at.next();
			for(int k=0;k<at.columns;++k)
			{
				col[k][i]=at.x[k];
			}
		}
	}
	~Table();
	void print()
	{
		cout<<"rows="<< col[0].size()<<" columns="<< col.size()<<endl;
		for(size_t k=0; k<col[0].size();++k)
		{
			cout<<"[ ";
			for(size_t j=0; j<col.size();++j)
				cout<<setw(8)<<col[j][k]<<" ";
			cout<<"]"<<endl;
		}
		cout<<endl;

	}
	vector<vector<double> > col;
};


class TableInterpolator: public GridTableInterpolatorBase
{
public:
	TableInterpolator(string dir):GridTableInterpolatorBase(dir)
	{
		Table temp;
		for(size_t i=0;i<file_list.size();++i)
		{
			string s=dir+file_list[i];
			tabs.push_back(temp);
			tabs[i].readFromFile(s.c_str());
		}
		value.resize(tabs[0].col.size()-1,0.0);
	}
	~TableInterpolator();
	void interpol(const double* Pos)
	{
//		print();
		setNeighbors(Pos);
		for(size_t ii=0;ii<value.size();++ii)
			value[ii]=0.0;

		for(size_t i=0;i<nv.size();++i)
			for(size_t ii=0;ii<value.size();++ii)
			{
//				cout<<i<<" "<<ii<<" "<<nv[i].prob<<" "<<interpolate(tabs[nv[i].ig].col[0],tabs[nv[i].ig].col[ii+1],Pos[keys.size()])<<endl;
				value[ii]+=nv[i].prob*interpolate(tabs[nv[i].ig].col[0],tabs[nv[i].ig].col[ii+1],Pos[keys.size()]);
			}

	}
	void print()
	{
		printBase();
		cout<<"Tables"<<endl;
		for(size_t i=0;i<tabs.size();++i)
			tabs[i].print();
	}
protected:
	vector<Table> tabs;
};


class TableInterpolatorBC: public TableInterpolator
{
public:
	TableInterpolatorBC(string dir):TableInterpolator(dir)
	{
		dir+="teff_list.dat";
		AsciiTable at(dir.c_str());
		teff.resize(at.rows,0.0);
		for(int i=0;i<at.rows;++i)
		{
			at.next();
			teff[i]=at.x[0];
		}

	};
	~TableInterpolatorBC();
	void interpol(const double* Pos)
	{
		setNeighbors(Pos);
		for(size_t ii=0;ii<value.size();++ii)
			value[ii]=0.0;

		int jj=locate(teff,Pos[keys.size()]);
		double temp1=1-(log10(Pos[keys.size()])-log10(teff[jj]))/(log10(teff[jj+1])-log10(teff[jj]));
//		cout<<teff[jj]<<" "<<Pos[keys.size()+1]<<endl;
		for(size_t i=0;i<nv.size();++i)
			for(size_t ii=0;ii<value.size();++ii)
			{
//				if(teff[jj]==32000.0)
//				{
//					cout<<i<<" "<<ii<<" "<<nv[i].prob<<" "<<nv[i].ig<<" "<<interpolate(tabs[nv[i].ig].col[0],tabs[nv[i].ig].col[ii+1],(teff[jj]+Pos[keys.size()+1]))<<endl;
//					cout<<i<<" "<<ii<<" "<<nv[i].prob<<" "<<nv[i].ig<<" "<<interpolate(tabs[nv[i].ig].col[0],tabs[nv[i].ig].col[ii+1],(teff[jj+1]+Pos[keys.size()+1]))<<endl;
//					cout<<(teff[jj]+Pos[keys.size()+1])<<endl;
//				}
				value[ii]+=nv[i].prob*interpolate(tabs[nv[i].ig].col[0],tabs[nv[i].ig].col[ii+1],teff[jj]+Pos[keys.size()+1])*temp1;
				value[ii]+=nv[i].prob*interpolate(tabs[nv[i].ig].col[0],tabs[nv[i].ig].col[ii+1],teff[jj+1]+Pos[keys.size()+1])*(1-temp1);
			}

	}
protected:
	vector<double> teff;
};



class Grid
{
public:
	Grid(){};
	void initialize(int rows, int columns)
	{
		if(data.size()>0)
		{
			data.clear();
		}
		data.resize(rows*columns,0.0);

	}
	void readFromFile(const char* fname)
	{
		AsciiTable at(fname);
		initialize(at.rows,at.columns);
		int ii=0;
		for(int i=0;i<at.rows;++i)
		{
			at.next();
			for(int k=0;k<at.columns;++k)
			{
				data[ii]=at.x[k]; ii++;
			}
		}
	}
	~Grid();
	void print()
	{
		for(size_t k=0; k<data.size();++k)
			cout<<data[k]<<endl;
		cout<<endl;
	}
	vector<double> data;
};


class GridInterpolator: public GridTableInterpolatorBase
{
public:
	GridInterpolator(string dir):GridTableInterpolatorBase(dir)
	{
		Grid temp;
		for(size_t i=0;i<file_list.size();++i)
		{
			string s=dir+file_list[i];
			grids.push_back(temp);
			grids[i].readFromFile(s.c_str());
		}

	}
	~GridInterpolator();
	void interpol(const double* Pos)
	{
		setNeighbors(Pos);
		for(size_t ii=0;ii<value.size();++ii)
			value[ii]=0.0;
		for(size_t i=0;i<nv.size();++i)
			for(size_t ii=0;ii<value.size();++ii)
				value[ii]+=nv[i].prob*grids[ii].data[nv[i].ig];

	}
	void print()
	{
		printBase();
		cout<<"Grids"<<endl;
		for(size_t i=0;i<grids.size();++i)
			grids[i].print();
	}
private:
	vector<Grid> grids;
};



//class Interpolator_d:
//{
//public:
//	Interpolator_d()
//	{
//
//	}
//	void readFromFile(string fname,string tagname)
//	{
//		string tname=tagname+".data";
//		fbvector<double> x(fname,"r",tname);
//		data.resize(x.size());
//		for(size_t i=0;i<x.size();++i)
//			data[i]=x[i];
//		ExtBinFormat2 ebf;
//		ebf.open(fname.c_str(),"r",tname.c_str(),0);
//		for(int8 i=0;i<ebf.rank();++i)
//
//
//	}
//	~GridInterpolator();
//	void interpol(const double* Pos)
//	{
//	}
//	void print()
//	{
//	}
//private:
//	vector<double> data;
//	vector<vector<double> > keys;
//
//};




#endif /* TABLEINTERPOLATOR_H_ */
