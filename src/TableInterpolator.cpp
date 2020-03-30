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

#include "TableInterpolator.h"

//TableInterpolator::TableInterpolator()
//{
//	// TODO Auto-generated constructor stub
//
//}

TableInterpolator::~TableInterpolator()
{
	// TODO Auto-generated destructor stub
}

TableInterpolatorBC::~TableInterpolatorBC()
{
	// TODO Auto-generated destructor stub
}


Table::~Table()
{
	// TODO Auto-generated destructor stub
}

Grid::~Grid()
{
	// TODO Auto-generated destructor stub
}


GridTableInterpolatorBase::~GridTableInterpolatorBase()
{

}

GridInterpolator::~GridInterpolator()
{
	// TODO Auto-generated destructor stub
}


AsciiTable::AsciiTable(const char* fname)
{
	_name="AsciiTable";
	char st1[256],st2[256];

	if((fd=fopen(fname,"r")))
	{

		if(check_comment())
			scan_status=fscanf(fd,"%s",st1);
		if(check_comment())
			scan_status=fscanf(fd,"%i%i%i%i",&rows,&dim1,&massid,&columns);
//			cout<<"rows="<<" "<<rows<<" columns="<<columns<<" ("<<dim1<<") massid="<<massid<<endl;
		if(columns>0)
		{
		dimid.resize(columns,0);
		for(int i=0;i<columns;++i)
		{
				scan_status=fscanf(fd,"%i",&dimid[i]);
			if(scan_status==0)
			{
				cout<<"Input format not correct, scan_status="<<scan_status<<endl;
				exit(1);
			}
//				cout<<i<<" "<<dimid[i]<<endl;
		}
		}
		else
		{
			columns=dim1;
			dimid.resize(columns);
			for(int i=0;i<columns;++i)
				dimid[i]=i;
		}

//		check duplicate dims
		{
			vector<int> dimid1=dimid;
			sort(dimid1.begin(),dimid1.end());
//				printv(dimid1.begin(),dimid1.end());
			for(size_t i=1;i<dimid1.size();++i)
				if(dimid1[i]==dimid1[i-1])
				{
					cout<<"Input format not correct duplicate dims"<<endl;
					exit(1);
				}

		}

		if(check_comment())
			scan_status=fscanf(fd,"%s",st2);
		if((strcmp(st1,st2)!=0)||(massid >= dim1)||((*max_element(dimid.begin(),dimid.end()))>=dim1)||(columns>dim1))
		{
			if((massid!=-1)&&((columns+1)>dim1))
			{
				cout<<"Input format not correct "<<massid<<" "<<columns<<" "<<dim1<<endl;
				exit(1);
			}

			cout<<"Input format not correct"<<endl;
			cout<<"<head>"<<endl;
			cout<<"rows cols massid dim dim_1 dim_2 ..dim_k"<<endl;
			cout<<"<head>"<<endl;
			cout<<"k should be equal to dim or dim should be 0"<<endl;
			cout<<"massid should be an integer or -1"<<endl;
			cout<<rows<<" "<<dim1<<" "<<massid<<" "<<columns<<" "<<(*max_element(dimid.begin(),dimid.end()))<<endl;
			exit(1);
		}
//			cout<<"Dims "<<columns<<" dim_id "<<endl;

		tempv.resize(dim1,0.0);
		x.resize(columns,0.0);
		ic=0;
	}
	else
	{
		cout<<"File not found "<<fname<<endl;
		exit(1);
	}
}
bool AsciiTable::check_comment()
{
	char *cptr;
//		for(int k=0;k<dim1;++k)
//			scan_status=fscanf(fd,"%f",&tempv[k]);
//		for(int k=0;k<columns;++k)
//			x[k]=tempv[dimid[k]];
	while(1)
	{
		fpos_t position;
		fgetpos (fd, &position);
		int c=fgetc(fd);
		fsetpos (fd, &position);
		if(feof(fd))
			return 0;
		if(c=='#')
			cptr=fgets(buf,1024,fd);
		else
			return 1;
	}
}
void AsciiTable::next()
{
	if(check_comment())
	{
		for(int k=0;k<dim1;++k)
			scan_status=fscanf(fd,"%f",&tempv[k]);
		for(int k=0;k<columns;++k)
			x[k]=tempv[dimid[k]];
	}
}
AsciiTable::~AsciiTable()
{
	if(fd!=NULL)
		fclose(fd);
}


GridTableInterpolatorBase::GridTableInterpolatorBase(string fname)
{
	fname+="interp_keys.txt";
	FILE* fd;
	char* cptr;
	char buf[1024];
	if((fd=fopen(fname.c_str(),"r")))
	{
		int i=0;
		while(1)
		{
			i++;
			cptr=fgets(buf,1024,fd);
			if(feof(fd))
				break;
			char* c1=NULL;
			c1=strchr(buf,'=');
			if(c1==NULL)
			{
				cout<<" records not found in"<<fname<<endl;
				exit(1);
			}

			if(strncmp(buf, "files=",c1-buf) == 0)
			{
				file_list.clear();
				stringSplit(c1+1," ,'\"\n",file_list);
			}
			else if(strncmp(buf, "keys=",c1-buf) == 0)
			{
				keys.push_back(vector<double>());
				vector<string> sv;
				stringSplit(c1+1," ,",sv);
				for(size_t k=0; k<sv.size();++k)
					keys.back().push_back(atof(sv[k].c_str()));
			}
		}
		fclose(fd);

	}
	else
	{
		cout<<"File not found "<<fname<<endl;
		exit(1);
	}

	dkeys.resize(keys.size(),0.0);
	for(size_t j=0; j<keys.size();++j)
	{
		dkeys[j]=keys[j][1]-keys[j][0];
		for(size_t k=1; k<keys[j].size();++k)
			if((keys[j][k]-keys[j][k-1]) != dkeys[j])
			{
				dkeys[j]=0.0;
				break;
			}
	}
	iv.resize(keys.size(),0);
	cPos.resize(keys.size(),0);
	du.resize(keys.size(),0);
	{
//			Neighbor ngb;
		nv.resize(int(pow(2.0,1.0*keys.size())),Neighbor());
	}

}
void GridTableInterpolatorBase::printBase()
{
	cout<<"file_list=[ ";
	for(size_t k=0; k<file_list.size();++k)
		cout<<file_list[k]<<" ";
	cout<<"]"<<endl;
	cout<<"Number of keys="<< keys.size()<<endl;
	for(size_t j=0; j<keys.size();++j)
	{
		cout<<"keys=[ ";
		for(size_t k=0; k<keys[j].size();++k)
			cout<<keys[j][k]<<" ";
		cout<<"]"<<endl;
	}
	cout<<endl;
}
void GridTableInterpolatorBase::setNeighbors(const double* Pos)
{
	for(size_t j=0;j<keys.size();++j)
	{
		cPos[j]=Pos[j];
		if(cPos[j]>=keys[j].back())
			cPos[j]=keys[j].back()*0.99999;
		if(cPos[j]<=keys[j].front())
			cPos[j]=keys[j].front()*1.00001;
		if (dkeys[j]>0.0)
			iv[j]=int((cPos[j]-keys[j][0])/(keys[j][1]-keys[j][0]));
		else
			iv[j]=locate(keys[j],Pos[j]);
		du[j]=(cPos[j]-keys[j][iv[j]])/(keys[j][iv[j]+1]-keys[j][iv[j]]);
	}


	double temp1,temp2,temp3;
	vector<int> ivc=iv;
	int ii=0;
	for(int i=iv[0];i<iv[0]+2;++i)
	{
		temp1=(i==iv[0])?(1-du[0]):du[0];
		ivc[0]=i;
		if(keys.size()>=2)
		{
			for(int j=iv[1];j<iv[1]+2;++j)
			{
				temp2=(j==iv[1])?(1-du[1]):du[1];
				ivc[1]=j;
				if(keys.size()==3)
				{
					for(int k=iv[2];k<iv[2]+1;++k)
					{
						temp3=(k==iv[2])?(1-du[2]):du[2];
						ivc[2]=k;
						nv[ii].prob=temp1*temp2*temp3;
						nv[ii].ig=index(&ivc[0]);
						ii++;
					}
				}
				else
				{
					nv[ii].prob=temp1*temp2;
					nv[ii].ig=index(&ivc[0]);
					ii++;
				}
			}
		}
		else
		{
			nv[ii].prob=temp1;
			nv[ii].ig=index(&ivc[0]);
			ii++;
		}
	}

}




