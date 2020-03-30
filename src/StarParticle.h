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


#ifndef STARPARTICLE_H_
#define STARPARTICLE_H_
#include "Functions.h"
#include<string>
#include<vector>



class StarParticlef
{
public:
	float Pos[17];
    inline float& feh() {return Pos[6];}
    inline float& alpha() {return Pos[7];}
	inline float& age() {return Pos[8];}
    inline float& smass() {return Pos[9];}
    inline float& rad() {return Pos[10];}
    inline float& mag(int k) {return Pos[11+k];}
   inline int& satID() {return *(int *)(Pos+14);}
   inline int& popID() {return *(int *)(Pos+15);}
   inline int& partID() {return *(int *)(Pos+16);}
};



class StarParticle
{
private:
	double Pos[21];
public:
    inline double& pos(int i) {return Pos[i];}
    inline double& feh() {return Pos[6];}
    inline double& alpha() {return Pos[7];}
	inline double& age() {return Pos[8];}
    inline double& smass() {return Pos[9];}
    inline double& rad() {return Pos[10];}
    inline double& mag(int k) {return Pos[11+k];}
    inline double& lum() {return Pos[14];}
    inline double& teff() {return Pos[15];}
    inline double& grav() {return Pos[16];}
    inline double& dmod() {return Pos[17];}
    inline int64_t& satID() {return *(int64_t *)(Pos+18);}
    inline int64_t& popID() {return *(int64_t *)(Pos+19);}
    inline int64_t& partID() {return *(int64_t *)(Pos+20);}
    void convert(StarParticlef &Starf)
    {
    	for(int i=0;i<14;++i)
    		Starf.Pos[i]=float(Pos[i]);
    	Starf.satID()=int(satID());
    	Starf.popID()=int(popID());
    	Starf.partID()=int(partID());
    }
    void print()
    {
    	cout<<setw(15)<<Pos[0]<<setw(15)<<Pos[3]<<setw(15)<<feh()<<setw(15)<<age()<<setw(15)<<smass()<<" "<<mag(0)<<" "<<mag(1)<<" "<<mag(2)<<" "<<dmod()<<endl;
    }
    inline void AbsToAppMag()
    {
    	rad()=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]+Pos[2]*Pos[2]);
    	dmod()=5*log10(rad()*100.0);
    	mag(0)+=dmod();
    	mag(1)+=dmod();
    	mag(2)+=dmod();
    }
    inline void AddPhotoError(Normaldev &gauss,int photoError,int magswap)
    {
    	double sigma_mag=0.0;
    	if(photoError==2)
    	{
    		double temp=pow(10.0,0.4*(min(25.0,mag(0))-24.5));
    		sigma_mag=sqrt((0.04-0.039)*temp+0.039*temp*temp);
    	}
    	else if(photoError==3)
    	{
    		double temp=pow(10.0,0.4*(min(28.0,mag(0))-27.5));
    		sigma_mag=sqrt((0.04-0.039)*temp+0.039*temp*temp);
    	}
    	else if(photoError==4)
    	{
    		double temp=pow(10.0,0.4*(min(23.0,mag(0))-22.6));
    		sigma_mag=sqrt((0.04-0.039)*temp+0.039*temp*temp);
    	}

		if(sigma_mag>0.0)
		{
			mag(0)+=gauss.dev()*sigma_mag;
			mag(1)+=gauss.dev()*sigma_mag;
			mag(2)+=gauss.dev()*sigma_mag;
			if(magswap==2)
				mag(2)=mag(0);
			else if (magswap==1)
				mag(1)=mag(0);
		}

    }

};


class ParticleTag
{
public:
	ParticleTag(string name1,int dim1,int id1,int datatype1,int status1):name(name1),dim(dim1),id(id1),status(status1), datatype(datatype1){}
	string name;
	int dim;
	int id;
	int status;
	int datatype;
	bool operator== (const ParticleTag& a)
	{
		return (a.name==name);
	}
};

class ParticleTag_eq: public unary_function<ParticleTag,bool>
{
	string s;
	int status;
public:
	explicit ParticleTag_eq(const string& ss,const int &st): s(ss), status(st){}
	bool operator() (const ParticleTag &t) const {return ((t.name==s)&&(t.status==status));}
};


class ParticleStar
{
public:
//	ParticleStar();
//	~ParticleStar();
    double* Pos;       // particle position
    inline double& FeH() {return Pos[6];}
    inline double& Alpha() {return Pos[7];}
    inline double& Age() {return Pos[8];}
    inline double& Mass() {return Pos[9];}
    inline double& h(int k) {return Pos[10+k];}
//    inline double& Density() {return Pos[11];}
     void getTags(vector<ParticleTag> &tags,int hdim1)
    {
    	int c=0;
    	int dtype=5;
    	tags.clear();
//    	tags.push_back(ParticleTag("/Pos",pos_size,c,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/Pos3",3,c,dtype,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/Vel3",3,c,dtype,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/FeH",1,c,dtype,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/Alpha",1,c,dtype,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/Age",1,c,dtype,1));c+=tags.back().dim;
    	tags.push_back(ParticleTag("/Mass",1,c,dtype,1));c+=tags.back().dim;
    	if(hdim1==3)
			{tags.push_back(ParticleTag("/H_smooth",2,c,dtype,1));c+=tags.back().dim;}
    	else if(hdim1==6)
    		{tags.push_back(ParticleTag("/H_cubic",2,c,dtype,1));c+=tags.back().dim;}
    	else
    	{cout<<"hdim must be 3 or 6"<<endl; exit(1);}
 //   	tags.push_back(ParticleTag("/Density",1,c,0));c+=tags.back().dim;
    }

};



//class Particle2 {
//public:
//	Particle2();
//	virtual ~Particle2();
//};

#endif /* PARTICLE2_H_ */
