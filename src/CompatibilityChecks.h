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

#ifndef COMPATIBILITYCHECKS_H_
#define COMPATIBILITYCHECKS_H_

#include "ebf.hpp"


void allCompatibilityChecks()
{
	char x1=0;
	int32_t x2=0;
	int64_t x3=0;
	float x4=0;
	double x5=0;
	int16_t x6=0;
	int8_t x9=0;
	uint8_t x10=0;
	uint16_t x11=0;
	uint32_t x12=0;
	uint64_t x13=0;
	int error=0;
	int status=1;
	if(sizeof(float)!=4)
		{cout<<"float not 32 bit"<<endl; status=0;}
	if(sizeof(int32_t)!=4)
		{cout<<"int not 32 bit"<<endl; status=0;}
	if(sizeof(double)!=8)
		{cout<<"double not 64 bit"<<endl; status=0;}
	if(sizeof(int64_t)!=8)
		{cout<<"int64 not 64 bit"<<endl; status=0;}
	if(sizeof(uint64_t)!=8)
		{cout<<"uint64 not 64 bit"<<endl; status=0;}
	if(status==0)
	{
		cout<<"EXIT(): 32 bit and 64 bit data type compatibility problem"<<endl;
		cout<<"hint: redefine datatypes"<<endl;
		error=1;
	}


	status=0;
	if((ebf::TypeV(x1)==1)&&(ebf::TypeV(x2)==2)&&(ebf::TypeV(x3)==3)&&(ebf::TypeV(x4)==4)&&(ebf::TypeV(x5)==5)&&(ebf::TypeV(x6)==6)&&(ebf::TypeV(x9)==9)&&(ebf::TypeV(x10)==10)&&(ebf::TypeV(x11)==11)&&(ebf::TypeV(x12)==12)&&(ebf::TypeV(x13)==13))
		status=1;
	if(status==0)
	{
		cout<<ebf::TypeV(x1)<<" "<<ebf::TypeV(x2)<<" "<<ebf::TypeV(x3)<<" "<<ebf::TypeV(x4)<<" "<<ebf::TypeV(x5)<<" "<<ebf::TypeV(x6)<<" "<<ebf::TypeV(x9)<<" "<<ebf::TypeV(x10)<<" "<<ebf::TypeV(x11)<<" "<<ebf::TypeV(x12)<<" "<<ebf::TypeV(x13)<<endl;
		cout<<"EXIT(): function getType not working properly"<<endl;
		error=1;
		exit(1);
	}
	if(error != 0)
		exit(1);

}



#endif /* COMPATIBILITYCHECKS_H_ */
