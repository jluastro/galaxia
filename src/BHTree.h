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

#ifndef BHTREE_H_
#define BHTREE_H_
//#include <functional>
#include "Population.h"
#include "ebf.hpp"



class BHNode
{
public:
	double Pos[4],l,mass,delta_mass,rho_max,rho_min,dmod,lage;
	BHNode* left;
	BHNode()
//	void init()
	{
		Pos[0]=0.0;
		Pos[1]=0.0;
		Pos[2]=0.0;
		Pos[3]=0.0;
		l=0.0;
		mass=0.0;
		delta_mass=0.0;
		rho_max=0.0;
		rho_min=0.0;
		dmod=0.0;
		lage=0.0;
		left=NULL;
	}
	~BHNode();
	void print()
	{
		cout<<"Pos=("<<Pos[0]<<","<<Pos[1]<<","<<Pos[2]<<") Mass="<<mass<<" delta_mass="<<delta_mass<<" rho_max="<<rho_max<<" rho_min="<<rho_min<<endl;
	}
	void scale_mass(double fac)
	{
		mass*=fac;
		delta_mass*=fac;
		rho_max*=fac;
		rho_min*=fac;
	}
//	inline bool operator<( const BHNode &a ) const    {	return dmod < a.dmod;    }
//	inline bool operator>( const BHNode &a ) const    {	return dmod > a.dmod;    }
	bool operator==( const BHNode &a ) const
	{
		if((Pos[0]!=a.Pos[0])||(Pos[1]!=a.Pos[1])||(Pos[2]!=a.Pos[2]))
			return 0;
		if((l!=a.l)||(mass!=a.mass)||(rho_max!=a.rho_max))
			return 0;
		if((rho_min!=a.rho_min)||(dmod!=a.dmod))
			return 0;
//		if((left==NULL)||(a.left==NULL))
		if(left!=a.left)
			return 0;
		return 1;
	}
	void pack(vector<double> &x,const BHNode* root)
	{
		if(x.size()!=12)
			x.resize(12);
		x[0]=Pos[0]; x[1]=Pos[1]; x[2]=Pos[2];x[3]=Pos[3];
		x[4]=l;
		x[5]=mass;
		x[6]=delta_mass;
		x[7]=rho_max;
		x[8]=rho_min;
		x[9]=dmod;
		x[10]=lage;
		*(int64_t* )&x[11]=left-root;
		if(left==NULL)
			*(int64_t* )&x[11]=-1;
	}
	void unpack(vector<double> &x,BHNode* root)
	{
		assert(x.size()==12);
		Pos[0]=x[0];
		Pos[1]=x[1];
		Pos[2]=x[2];
		Pos[3]=x[3];
		l=x[4];
		mass=x[5];
		delta_mass=x[6];
		rho_max=x[7];
		rho_min=x[8];
		dmod=x[9];
		lage=x[10];
		left=root+(*(int64_t* )&x[11]);
		if((*(int64_t* )&x[11])==(-1))
			left=0;
	}
	bool splittable(const double &mass_split,const double &l_split)
	{
//		return (((mass > mass_split) && (l > l_split))||(l > 5.0) );
		return ((mass > mass_split)||(l > 5.0));
	}
	void setUpRoot(const double* Pos1,double l1,Population &Pop)
	{
		left=NULL;
		l=l1;
		lage=Pop.dage;
		if(Pop.optionE==0)
			lage=0.0;
	   	Pos[0]=Pos1[0];
	   	Pos[1]=Pos1[1];
	   	Pos[2]=Pos1[2];
	   	Pos[3]=Pop.age;
	   	computeMass(Pop);
	   	cout<<"Age="<<Pos[3]<<" lAge="<<lage<<endl;
	}
	void initializeFromParent(const double* Pos1,double l1,double lage1,int i1,int i2,int i3,Population &Pop)
	{
		left=NULL;
		l=l1/2;
		lage=lage1;
	   	Pos[0]=Pos1[0]+l*i1;
	   	Pos[1]=Pos1[1]+l*i2;
	   	Pos[2]=Pos1[2]+l*i3;
	   	Pos[3]=Pos1[3];
	   	computeMass(Pop);
	}
	void computeMass(Population &Pop)
	{
//		double xmind[3],xmaxd[3];
		Pop.computeMinMaxDensity(Pos,l,lage,rho_min,rho_max);
//		rho_min=Pop.density(xmind);
//		rho_max=Pop.density(xmaxd);
		delta_mass=fabs(1-rho_min/rho_max);
		if(lage>0)
			mass=rho_max*l*l*l*lage*16;
		else
			mass=rho_max*l*l*l*8;

//		cout<<xmax[0]<<" "<<xmax[1]<<" "<<xmax[2]<<endl;
//		cout<<mass<<" "<<delta_mass<<" "<<l<<endl;
	}
	void calculateDmod(const double* Pos1)
	{
		dmod=(Pos[0]-Pos1[0])*(Pos[0]-Pos1[0])+(Pos[1]-Pos1[1])*(Pos[1]-Pos1[1])+(Pos[2]-Pos1[2])*(Pos[2]-Pos1[2]);
		dmod=5*log10(sqrt(dmod)*100.0);
	}
	void generatePos(StarParticle &Star,Population &Pop)
	{
		double temp=0,tempr=1;
		while(tempr>temp)
		{
			Star.pos(0)=Pos[0]+(grandomu()-0.5)*2*l;
			Star.pos(1)=Pos[1]+(grandomu()-0.5)*2*l;
			Star.pos(2)=Pos[2]+(grandomu()-0.5)*2*l;
			if(Pop.optionE==1)
				Star.age()=Pos[3]+(grandomu()-0.5)*2*lage;
			temp=Pop.density(&Star.pos(0),Star.age());
			tempr=grandomu()*rho_max*2.0;
		}
//		Pop.generateVel(Pos1);
	}
	void computeMassRefined(double tol,int maxref,Population &Pop)
	{
		double x[4],temp,h,temp1,temp2;//,temp3;
		int n=1,ref;
		vector<double> mass1(maxref+1);
		mass1[0]=mass;

		temp1=Pop.density(Pos,Pos[3]);
		temp2=Pop.density(Pos,Pos[3]);

//		int option=1;
//		if(option==1)
		for(ref=1;ref<=maxref;++ref)
		{
			temp=0.0;
			h=l/n;
			n*=2;
			for(int i=0;i<n;++i)
			{
				  x[0]=(Pos[0]-l)+(i+0.5)*h;
//				  cout<<"x0 "<<x[0]<<endl;
				for(int j=0;j<n;++j)
				{
					  x[1]=(Pos[1]-l)+(j+0.5)*h;
//					  cout<<"x1 "<<x[1]<<endl;
					for(int k=0;k<n;++k)
					{
						  x[2]=(Pos[2]-l)+(k+0.5)*h;
//						  cout<<"x2 "<<x[2]<<endl;

						  if(lage>0.0)
						  {
							for(int l1=0;l1<n;++l1)
							{
								  x[3]=(Pos[3]-lage)+(l1+0.5)*lage*2/n;
								  temp+=Pop.density(x,x[3]);
//								  temp3=Pop.density(x,x[3]);
							}
						  }
						  else
						  {
							  temp+=Pop.density(x,x[3]);
							  //							  temp3=Pop.density(x,x[3]);
						  }

//						  if(temp2<temp3)
//							  temp2=temp3;
//						  else if(temp1>temp3)
//							  temp1=temp3;
					}
				}
			}

			if(lage>0.0)
				mass1[ref]=temp*h*h*h*lage*2/n;
			else
				mass1[ref]=temp*h*h*h;


//			cout<<"Ref "<<ref<<" "<<mass<<" "<<mass1[ref]<<" "<<(mass1[ref-1]-mass1[ref])/mass1[ref]<<endl;
			if((mass1[ref]<1e-30)&&(ref>1))
				break;
			if((fabs(mass1[ref-1]-mass1[ref]))<(0.5*sqrt(0.5*(mass1[ref]+mass1[ref-1]))) )
//			if((fabs(mass1[ref-1]-mass1[ref]))<tol)
				break;
		}

		if(ref<maxref)
		{
			delta_mass=fabs(mass1[ref]-mass1[ref-1])/mass1[ref];
			mass=mass1[ref];
		}
		else
		{
			delta_mass=fabs(mass1[maxref]-mass1[maxref-1])/mass1[maxref];
			mass=mass1[maxref];
		}

//		temp1=computeMassRec(Pos,l,lage,Pop);
//		mass=temp1;

//		if(fabs(temp1-mass)/sqrt(temp1)>0.5)
//			cout<<mass<<" "<<temp1<<" "<<ref<<" "<<fabs(temp1-mass)/sqrt(temp1)<<endl;

//		cout<<"Ref "<<ref<<" "<<fabs(mass1[ref-1]-mass1[ref])<<" "<<fabs(mass1[ref-1]-mass1[ref])*2/(mass1[ref]+mass1[ref-1])<<" "<<(mass1[ref]+mass1[ref-1])<<endl;
//		if((0.5*(mass1[ref]+mass1[ref-1]))>100.0)
//		cout<<"Ref "<<ref<<" Mass"<<(0.5*(mass1[ref]+mass1[ref-1]))<<" deltaM/sqrt(M)"<<fabs(mass1[ref-1]-mass1[ref])*2/sqrt(0.5*(mass1[ref]+mass1[ref-1]))<<" "<<rho_max/temp2<<endl;

	}

	double computeMassRec(double* Pos1,double l1,double lage1,Population &Pop)
	{
		double x[4],temp,temp1,temp2;
		temp=0.0;
		double vol=l1*l1*l1*lage1;
		temp2=Pop.density(Pos1,Pos1[3])*vol;
		for(int i=0;i<2;++i)
		{
			x[0]=(Pos1[0]-l1)+(i+0.5)*l1;
			for(int j=0;j<2;++j)
			{
				x[1]=(Pos1[1]-l1)+(j+0.5)*l1;
					for(int k=0;k<2;++k)
					{
						x[2]=(Pos1[2]-l1)+(k+0.5)*l1;
						if(lage1>0.0)
						{
							for(int p=0;p<2;++p)
							{
								x[3]=(Pos1[3]-lage1)+(p+0.5)*lage1;
								temp1=Pop.density(x,x[3])*vol;
								if(fabs(temp1-temp2)>0.5*sqrt(0.5*(temp1+temp2)))
									temp+=computeMassRec(x,l1/2,lage1/2,Pop);
								else
									temp+=temp1;
							}
						}
						else
						{
							temp1=Pop.density(x,x[3])*vol;
							if((temp1)>1.0)
								temp+=computeMassRec(x,l1/2,lage1/2,Pop);
							else
								temp+=temp1;
						}
					}
			}
		}

			if(lage1>0.0)
				return temp;
			else
				return temp;

	}
	inline BHNode& child(int k){return left[k];}
	void split(BHNode* Node1,Population &Pop)
	{
		left=Node1+1;
		left[0].initializeFromParent(Pos,l,lage,-1,-1,-1,Pop);
		left[1].initializeFromParent(Pos,l,lage,-1,-1, 1,Pop);
		left[2].initializeFromParent(Pos,l,lage,-1, 1,-1,Pop);
		left[3].initializeFromParent(Pos,l,lage,-1, 1, 1,Pop);
		left[4].initializeFromParent(Pos,l,lage, 1,-1,-1,Pop);
		left[5].initializeFromParent(Pos,l,lage, 1,-1, 1,Pop);
		left[6].initializeFromParent(Pos,l,lage, 1, 1,-1,Pop);
		left[7].initializeFromParent(Pos,l,lage, 1, 1, 1,Pop);
	}
};

class nodep_less_than: binary_function<BHNode*, BHNode*,bool>
{
public:
	bool operator() (BHNode* const& N1,BHNode* const& N2) const {return N1->dmod < N2->dmod;}
//	bool operator() (BHNode* const& N1,BHNode* const& N2) const {return N1->mass < N2->mass;}
};


class BHTree
{
public:
	BHTree()
	{
		outputFile="bhtree_";
	}
	BHTree(Population &Pop,const double* posC,int warpFlareOn1,int option,const string inputDir);
	void initialize(Population &Pop,const double* posC,int warpFlareOn1,int option,const string inputDir);
	void save(const string fname)
	{
		cout<<left<<setw(36)<<"Writing tree to file: "<<fname<<"....";
		vector<double> x(12,0.0);
		int64_t dims[2];
		dims[0]=rootv.size();
		dims[1]=x.size();
		ebf1.Open(fname,"/Nodes","w",5,"",2,dims);
		for(size_t i=0;i<rootv.size();++i)
		{
			rootv[i].pack(x,root);
//			fwrite(&x[0], sizeof(x[0]) * x.size(), 1, ebf1.fp());
			ebf1.Write(&x[0],x.size());
		}
		ebf1.Close();

		vector<int64_t> leafs_ind(leafs.size(),0);
		for(size_t i=0;i<leafs.size();++i)
			leafs_ind[i]=leafs[i]-root;
		ebf::Write(fname,"/Leafs",&leafs_ind[0],"a","",leafs.size());
		cout<<"Done"<<endl;

	}
	void load(const string fname)
	{

		cout<<left<<setw(36)<<"Reading tree from file- "<<fname<<endl;
//		cout<<setw(36)<<" "<<"Checking tree integrity--"<<flush;
//		ebfformat2::ebfcheck(fname);
//		cout<<"OK"<<endl;

		vector<double> x(12,0.0);
		ebf1.Open(fname,"/Nodes");
		assert(ebf1.ecode==0);
		assert(ebf1.dim(0)>0);
		assert(ebf1.dim(1)==12);
		rootv.resize(ebf1.dim(0));
		root=&rootv[0];
		for(size_t i=0;i<rootv.size();++i)
		{
			ebf1.Read(&x[0],x.size());
//			read(&x[0], sizeof(x[0]) * x.size(), 1, ebf1.fp());
			rootv[i].unpack(x,root);
		}
		ebf1.Close();


		ebf1.Open(fname,"/Leafs");
		assert(ebf1.ecode==0);
		assert(ebf1.dim(0)>0);
		vector<int64_t> leafs_ind;
		ebf::Read(fname,"/Leafs",leafs_ind);
		leafs.resize(leafs_ind.size());
		for(size_t i=0;i<leafs.size();++i)
			leafs[i]=root+leafs_ind[i];
	}
	~BHTree();
	ebf::EbfFile ebf1;
	vector<BHNode> rootv;
	BHNode* root;
//	char outputFile[500];
	string outputFile;
	vector<BHNode*> leafs;
};

#endif /* BHTREE_H_ */
