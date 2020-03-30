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

#include <utility>
#include <iostream>
#include "Satellite.h"
#include "Functions.h"
#include "ebfvector.hpp"

#define ALL(A) A.begin,A.end()
using namespace std;

vector<int> lindgen(int N)
{
	vector<int> x(N);
	for(int i=0;i<N;++i)
		x[i]=i;
	return x;
}



void Satellite::print()
{
	cout << "Satellite Info" << endl;
//	cout << "M_tot " << m_tot << endl;
//	cout << "Feh min max array size" << (*min_element(feh.begin(), feh.end()))
//			<< " " << (*max_element(feh.begin(), feh.end())) << " "
//			<< feh.size() << endl;
//	cout << "Feh min max array size" << (*min_element(feh_a.begin(),
//			feh_a.end())) << " " << (*max_element(feh_a.begin(), feh_a.end()))
//			<< " " << feh_a.size() << endl;
//	cout << "Age min max array size" << (*min_element(age.begin(), age.end()))
//			<< " " << (*max_element(age.begin(), age.end())) << " "
//			<< age.size() << endl;
//	cout << "Age min max array size" << (*min_element(age_a.begin(),
//			age_a.end())) << " " << (*max_element(age_a.begin(), age_a.end()))
//			<< " " << age_a.size() << endl;
//	cout<<"hdim"<<hdim<<endl;
//	printv(age_a.begin(), age_a.end(), "age_a");
//	printv(faca_a.begin(), faca_a.end(), "faca_a");
//	printv(feh_a.begin(), feh_a.end(), "feh_a");
	//	printv(m1_a.begin(), m1_a.end(), "m1_a");
	//	printv(m2_a.begin(), m2_a.end(), "m2_a");
	//	printv(fac_a.begin(), fac_a.end(), "fac_a");
	//				PSplot plot1;
	//				plot1.plot(age_sampler.x,age_sampler.cpd,"x","P(<x)");
}

void Satellite::initialize_stellar_data(Sampler& imf,IsochroneDB &ic)
{
	if(debug) cout <<"Satellite Initializing ....... ";
	age.resize(nsize);
	feh.resize(nsize);
	alpha.resize(nsize);
	int k = 0;
	for (ParticleStar* Part = pBegin; Part < pEnd; Part++)
	{
		age[k] = log10((Part->Age()) * 1.0e9);
//		feh[k] = 0.019 * pow(10.0, Part->FeH());
		feh[k] = Part->FeH();
		alpha[k] = Part->Alpha();
		m_tot += Part->Mass();
		k++;
	}
	m_tot *= 2.0;

	alpha2 = alpha;
	feh2 = feh;
	sort2(feh2.begin(), feh2.end(), alpha2.begin());
	vector<double>::iterator it1 = unique2(feh2.begin(), feh2.end(), alpha2.begin());
	alpha2.resize(it1 - feh2.begin());
	feh2.resize(it1 - feh2.begin());

	vector<int> ind=lindgen(age.size());
	sorti2(age.begin(), age.end(), ind.begin());
	orderv(feh,ind);
	orderv(alpha,ind);
	vector<double> age1 = age;
	vector<double> feh1 = feh;
	ageSampler = new Sampler(age);



	ind=lindgen(age.size());
	vector<double>::iterator it = unique2(age.begin(), age.end(), ind.begin());
	orderv(feh,ind);
	orderv(alpha,ind);
	feh.resize(it - age.begin());
	alpha.resize(it - age.begin());
	age.resize(it - age.begin());

	ind=lindgen(age.size());
	it = unique2(feh.begin(), feh.end(), ind.begin());
	orderv(age,ind);
	orderv(alpha,ind);
	age.resize(it - feh.begin());
	alpha.resize(it - feh.begin());
	feh.resize(it - feh.begin());


	linspace(age_a, (*min_element(age.begin(), age.end())), (*max_element(
			age.begin(), age.end())), 0.02);
//	cout<<(*min_element(age.begin(), age.end()))<<(*max_element(age.begin(), age.end()))<<endl;
	//	linspace(age_a,(*min_element(ALL(age))),(*max_element(ALL(age))),0.5);
	//	linspace(age_a,pow10.0,age_min-9),pow(10.0,age_max-9),0.5);
	//	transform(age_a.begin(),age_a.end(),log10);
	interpolate(age, feh, age_a, feh_a);
	//	ic.min_max_m(age_a,feh_a,feh_a,m1_a,m2_a,5.0);
	//	imf.getFacv(m1_a,m2_a,fac_a);
	interpolate(ageSampler->x, ageSampler->cpd, age_a, faca_a);


//	iso.resize(age1.size());
//	for (size_t i = 0; i < age1.size(); ++i)
//	{
//		iso[i] = ic.nearest_ic(age1[i], feh1[i], 0.0);
//	}
//	int count = 1;
//	iso1.push_back(iso[0]);
//	for (size_t i = 1; i < age1.size(); ++i)
//	{
//		if (iso[i - 1] == iso[i])
//		{
//			count++;
//		}
//		else
//		{
//			iso1_pro1.push_back(count * 1.0 / age1.size());
//			count = 1;
//			iso1.push_back(iso[i]);
//		}
//	}
//	iso1_pro1.push_back(count * 1.0 / age1.size());




	if(debug) cout << "Done" << endl;
}

Satellite::~Satellite()
{
	delete ageSampler;
	ageSampler = NULL;
	delete[] data1;
	delete[] part;
}

double kernel_func(double x)
{
	double x2=x*x;
	return (1-x2)*x2*x2*x;
}

void Satellite::spawn1(SurveyDesign &sur, Sampler& imf,IsochroneDB &ic,double fSample,int seed1)
{
	print();
	double distMod;
	initialize_stellar_data(imf, ic);
	Sampler kerSampler(1000,0.0,1.0,kernel_func,0);
	double pr[6];
	Normaldev randomn(0.0,1.0,seed1);
	Ran randomu(seed1+12);

	StarParticle Star,Star1;
	vector<double> m3_a;
	m3_a=m2_a;
	nstars = 0;
	nstars1 = 0;
	int nstars2 = 0;
	int nstars3 = 0;
	l_tot_v = 0.0;
	cout<<"Particles="<< pEnd - pBegin<<" Mass="<<m_tot <<" "<<imf.meanx<< endl;
	for (ParticleStar* Part = pBegin; Part < pEnd; Part++)
	{

		for(int i=0;i<6;++i)
			Star.pos(i)=Part->Pos[i]-sur.posC[i];
		distMod = distance(Part->Pos, sur.posC, 3) - Part->h(0);
		if((sur.geo->checkSphere(&(Star.pos(0)),Part->h(0)) == 1)&&(distMod<=sur.r_max))
		{
			nstars2++;
			distMod = 5 * log10(max(distMod,0.01) * 100);
			ic.min_max_m(age_a, feh_a, feh_a, m1_a, m2_a, min(sur.absMag[1],sur.appMag[1]- distMod)+0.6,sur.All->starType);
			imf.getFacv(m1_a, m2_a, fac_a);
//			printv(age_a.begin(), age_a.end(), "age_a");
//			printv(faca_a.begin(), faca_a.end(), "faca_a");
//			printv(m1_a.begin(), m1_a.end(), "m1_a");
//			printv(m2_a.begin(), m2_a.end(), "m2_a");
//			printv(fac_a.begin(), fac_a.end(), "fac_a");
//			exit(1);
			int status=0;
			for (unsigned int j = 1; j < age_a.size(); ++j)
			{

				double temp = Part->Mass() * 2 * fSample * fac_a[j] * (faca_a[j]- faca_a[j - 1]) / imf.meanx;

				//			print();
				int stars = int(temp);
				if (randomu.doub() <= (temp - stars))
					stars++;

				if (stars > 0)
				{
					ageSampler->setRange(age_a[j - 1], age_a[j]);
					imf.setRange(m1_a[j], m2_a[j]);
				}
				nstars3+=stars;
				for (int l = 0; l < stars; ++l)
				{
					for(int i=0;i<6;++i)
						Star.pos(i)=Part->Pos[i]-sur.posC[i];


					double temp=0.0,temp1=0.0;
					temp1=0.0;
					for(int i=0;i<hdim;++i)
					{
						pr[i]=randomn.dev();
						temp1+=pr[i]*pr[i];
					}
					temp1=sqrt(temp1);
					temp=kerSampler.rand();

					for(int i=0;i<3;++i)
						Star.pos(i)+=Part->h(0)*pr[i]*temp/temp1;
					if(hdim==6)
					{
						for(int i=0;i<3;++i)
							Star.pos(i+3)+=Part->h(1)*pr[i+3]*temp/temp1;
					}

					Star.age() = ageSampler->rand();
					interpolate(age, feh, Star.age(), Star.feh());
					Star.smass() = imf.rand();
					Star.alpha() = 0.0;
					Star.mag(0) = 0.0;
					Star.mag(1) = 0.0;
					Star.mag(2) = 0.0;
					Star.satID() = satID;
					Star.popID() = 10;
					Star.partID() = 1;

					ic.interpolateStar(Star);


					if(status==0)
					{
						Star1=Star;
						for(int i=0;i<6;++i)
							Star1.pos(i)=Part->Pos[i]-sur.posC[i];
						Star1.AbsToAppMag();
						Star1.AddPhotoError(randomn,sur.All->photoError,ic.magswap);
						interpolate(feh2, alpha2, Star1.feh(), Star1.alpha());
						Star1.partID() = 0;
					}

					Star.AbsToAppMag();
					Star.AddPhotoError(randomn,sur.All->photoError,ic.magswap);
					interpolate(feh2, alpha2, Star.feh(), Star.alpha());


//					cout<<nstars+nstars1<<" "<<nstars<<" "<<temp<<" "<<status<<endl;
//					Star.print();
//					Star1.print();
//					if((nstars+nstars1)==10)
//						exit(1);

// for testing scattering errors
//					Star.alpha()=1;
//					Star.pos(4)=Part->h(0)*temp;
//					Star.pos(5)=Part->h(1)*temp;
//					Star1.alpha()=0;
//					Star1.pos(4)=0.0;
//					Star1.pos(5)=0.0;

					if(status==0)
					{
						if (sur.push_check(Star))
						{
							if(sur.push_back(Star1))
							{
//								l_tot_v += pow(10.0, (4.83 - Star1.mag(0) + Star1.dmod())/ 2.5);
								nstars++;
								status=1;
							}
							else if(sur.push_back(Star))
							{
//								l_tot_v += pow(10.0, (4.83 - Star.mag(0) + Star.dmod())/ 2.5);
								nstars++;
							}
							else
								nstars1++;
						}
						else
							nstars1++;

					}
					else
					{
						if (sur.push_back(Star))
						{
//							l_tot_v += pow(10.0, (4.83 - Star.mag(0) + Star.dmod())/ 2.5);
							nstars++;
						}
						else
							nstars1++;
					}


					if (((nstars + nstars1) % 10000000) == 0)
					{
						cout<< nstars<< " accepted  " << nstars1<<" rejected Parts="<<nstars2<<" outof "<<Part-pBegin<<endl;
					}
				}
			}
		}
	}

	cout << "Total Stars="<<nstars3<<" accepted="<< nstars << " rejected=" << nstars1 << endl;
	sur.flush();
}



void Satellite::initialize(int size1)
{
	nsize = size1;
//	int	dim = 0;
	items = 0;
	for (size_t i = 0; i < tags.size(); ++i)
	{
//		if (tags[i].status)
//			dim += tags[i].dim;
		items += tags[i].dim;
	}
	//	items=items1;
	//	dim=dim1;
	data1 = new double[nsize * items];
	data = data1;
	part = new ParticleStar[nsize];
	for (int i = 0; i < nsize; i++)
		part[i].Pos = &data[i * items];

	vector<ParticleTag>::iterator it;
	it = find_if(tags.begin(), tags.end(), ParticleTag_eq("/H_smooth", 0));
	if (it != tags.end())
		for (long i = 0; i < long(nsize); i++)
		{
			part[i].Pos[it->id] = 1.0;
			part[i].Pos[it->id + 1] = 1.0;
		}
	pBegin = part;
	pEnd = part + nsize;
	if(debug) cout<<"Particles="<<nsize<<endl;
}

Satellite::Satellite(const string& fname, int sat_no, int satID1,int hdim1,int nres1)
{
	hdim=hdim1;
	nres=nres1;
	satID = satID1;
	pBegin = NULL;
	pEnd = NULL;
	items = 0;
	nsize = 0;
	m_tot = 0.0;
	nstars = 0;
	debug=1;
	//	initialize();
	ageSampler = NULL;
	ParticleStar s;
	s.getTags(tags, hdim);
	readEbfFile(fname, sat_no);
}




void Satellite::readEbfFile(const string &fname, int sat_no)
{
	ebf::EbfVector<float> fb_pos3(fname,"/Pos3");
	ebf::EbfVector<float> fb_vel3(fname,"/Vel3");
	ebf::EbfVector<float> fb_feh(fname,"/FeH");
	ebf::EbfVector<float> fb_alpha(fname,"/Alpha");
	ebf::EbfVector<float> fb_age(fname,"/Age");
	ebf::EbfVector<float> fb_mass(fname,"/Mass");
	initialize(fb_mass.size());

	{
		ebf::EbfVector<int> fb_satid(fname,"/id");
		satID=fb_satid[0];
	}



	for(int i=0;i<nsize;++i)
	{
		part[i].Pos[0]=fb_pos3[i*3+0];
		part[i].Pos[1]=fb_pos3[i*3+1];
		part[i].Pos[2]=fb_pos3[i*3+2];
		part[i].Pos[3]=fb_vel3[i*3+0];
		part[i].Pos[4]=fb_vel3[i*3+1];
		part[i].Pos[5]=fb_vel3[i*3+2];
		part[i].FeH()=fb_feh[i];
		part[i].Alpha()=fb_alpha[i];
		part[i].Mass()=fb_mass[i];
		part[i].Age()=fb_age[i];
	}


	string s(fname);
	if(hdim==3)
	{
		if(nres==32)
			s.insert(s.find(".ebf"),"_d3n32_den");
		else if(nres==64)
			s.insert(s.find(".ebf"),"_d3n64_den");
		else if(nres==128)
			s.insert(s.find(".ebf"),"_d3n128_den");
		else
			s.insert(s.find(".ebf"),"_d3n64_den");
	}
	else
		s.insert(s.find(".ebf"),"_d6n64_den");

	ebf::EbfVector<float> fb_hcubic(s,"/H_cubic");
	if(fb_hcubic.rank()==1)
	{
		for(int i=0;i<nsize;++i)
		{
			part[i].h(0)=fb_hcubic[i];
			part[i].h(1)=0.0;
		}
	}
	if(fb_hcubic.rank()==2)
	{
		for(int i=0;i<nsize;++i)
		{
			part[i].h(0)=fb_hcubic[i*2+0];
			part[i].h(1)=fb_hcubic[i*2+1];
		}
	}

}

