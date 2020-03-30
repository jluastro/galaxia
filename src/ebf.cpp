/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Author: Sanjib Sharma                                                     *
 * Copyright (c) 2012 Sanjib Sharma                                          *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of EBF.  The full EBF copyright notice, including       *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree of EBF.                              *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "ebf.hpp"
#include <cstring>
#include<cstdlib>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>


//---------------------------------------------------------------------
namespace ebf
{


template<> int TypeT<char>(void)                {return 1;}
template<> int TypeT<int32_t> (void)                {return 2;}
template<> int TypeT<int64_t> (void)               {return 3;}
template<> int TypeT<efloat32> (void)              {return 4;}
template<> int TypeT<efloat64> (void)             {return 5;}
template<> int TypeT<int16_t> (void)          {return 6;}
template<> int TypeT<std::string> (void)        {return 7;}
template<> int TypeT<int8_t> (void)        {return 9;}
template<> int TypeT<uint8_t> (void)      {return 10;}
template<> int TypeT<uint16_t> (void) {return 11;}
template<> int TypeT<uint32_t> (void)       {return 12;}
template<> int TypeT<uint64_t> (void)              {return 13;}

void Info(const std::string& filename1)
{
	using namespace std;
	ebfc::EbfHeader ebfh;
	cout << left << setw(28) << "name" << setw(9) << "dtype"<< setw(12) << "units" << left<< " dim" << endl;
	cout<< "--------------------------------------------------------------------------"<< endl;
	FILE* fd = fopen(filename1.c_str(), "rb");
	if (fd != NULL)
	{
		fseek(fd, 0, SEEK_END);
		int64_t loc = ftell(fd);
		fseek(fd, 0, SEEK_SET);
		while (ftell(fd) < loc)
		{
			ebfc::EbfHeader_read(&ebfh,fd);
			fseek(fd,ebfh.capacity_,SEEK_CUR);
			string s="("+string(ebfh.unit)+")";
			cout << setw(28) << left << ebfh.name << setw(9)
								<< ebfutils::EbfCodeToName(ebfh.datatype)<<setw(12)<<s<< left;

			cout << " (";
			for (int i = 0; i < ebfh.rank - 1; ++i)
				cout << ebfh.dim[i] << ",";
			cout << ebfh.dim[ebfh.rank - 1] << ")" <<" "<<ftell(fd)<< endl;
			if(ebfh.datatype==8)
			{
				cout<<"structure definition:"<<endl;
				cout<<ebfh.sdef<<endl;
			}
			if(ebfh.ecode!=0)
			{
				cout<<"Ebf error: in reading header"<<endl;
				exit(1);
			}
		}
		fclose(fd);
	}
	else
	{
		cout<<"Ebf error: file not found"<<endl;
	}
}


int Copy(const std::string& infile, const std::string& outfile,const std::string& mode, const std::string& dataname)
{
	return ebfc::Ebf_Copy(infile.c_str(),outfile.c_str(),mode.c_str(),dataname.c_str());
}

//--------------------------------------------------------------------------------
namespace ebfutils
{
void Ebf_Error(const std::string& mystr)
{
	throw std::logic_error(mystr);
}

int64_t FileSize(const std::string& filename)
{
	int64_t size1=0;
	FILE* fp = fopen(filename.c_str(), "r");
	if (fp != NULL)
	{
		fseek(fp, 0, SEEK_END);
		size1 = ftell(fp);
		rewind(fp);
		fclose(fp);
	}
	return size1;
}
bool FileExists(const std::string& filename)
{
	FILE* fp=fopen(filename.c_str(),"r");
	if(fp!=NULL)
	{
		fclose(fp);
		return 1;
	}
	else
	{
		return 0;
	}
}

void SwapEndian(void* addr, int datasize,int64_t nmemb)
{
	int64_t pattern[]={datasize,nmemb,0};
	int64_t i, j = 0, k;
	char c;
	while (pattern[j] > 0)
	{
		for (k = 0; k < pattern[j + 1]; k++)
		{
			for (i = 0; i < pattern[j] / 2; i++)
			{
				c = *((char*) addr + i);
				*((char*) addr + i) = *((char*) addr + (pattern[j] - i - 1));
				*((char*) addr + (pattern[j] - i - 1)) = c;
			}
			addr = ((char *) addr) + pattern[j];
		}
		j += 2;
	}
}
int TypeSize(int x)
{
	switch (x)
	{
	  case 1:
		  return 1;
	    break;
	  case 2:
		  return 4;
	    break;
	  case 3:
		  return 8;
	    break;
	  case 4:
		  return 4;
	    break;
	  case 5:
		  return 8;
	    break;
	  case 6:
		  return 2;
	    break;
	  case 9:
		  return 1;
	    break;
	  case 10:
		  return 1;
	    break;
	  case 11:
		  return 2;
	    break;
	  case 12:
		  return 4;
	    break;
	  case 13:
		  return 8;
	    break;
	  default:
		  throw std::runtime_error("Ebf Error from ebfutils::getDataSize()-unrecognized data type ");
//		  exit(1);
	    return 0;
	  }
}






std::string EbfCodeToName(int x)
{
	std::string s;

	switch (x)
	{
	  case 1:
		  s="char";
	    break;
	  case 2:
		  s="int32";
	    break;
	  case 3:
		  s="int64";
	    break;
	  case 4:
		  s="float32";
	    break;
	  case 5:
		  s="float64";
	    break;
	  case 6:
		  s="int16";
	    break;
	  case 7:
		  s="string";
	    break;
	  case 8:
		  s="struct";
	    break;
	  case 9:
		  s="int8";
	    break;
	  case 10:
		  s="uint8";
	    break;
	  case 11:
		  s="unit16";
	    break;
	  case 12:
		  s="uint32";
	    break;
	  case 13:
		  s="uint64";
	    break;
	  default:
		  s="undefined";
		  break;
	  }
	return s;
}


void StringSplit(const std::string &s,const char* delimiters,std::vector<std::string> &sv)
{
	sv.clear();
	int length = s.size();
	char* c1= new char[length+1];
	strncpy(c1,s.data(),length);
	c1[length]=0;

	char * pch;
	pch = strtok (c1,delimiters);
	while (pch != NULL)
	{
		sv.push_back(pch);
	    pch = strtok (NULL,delimiters);
	}

	delete[] c1;
}

}


//------------------------------------------------------

static void ebf_tfunc_char_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}
static void ebf_tfunc_char_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(char *)(data1+i*stride1);
}


static void ebf_tfunc_int32_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}
static void ebf_tfunc_int32_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(int32_t *)(data1+i*stride1);
}


static void ebf_tfunc_int64_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}
static void ebf_tfunc_int64_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(int64_t *)(data1+i*stride1);
}


static void ebf_tfunc_float32_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}
static void ebf_tfunc_float32_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(float *)(data1+i*stride1);
}


static void ebf_tfunc_float64_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}
static void ebf_tfunc_float64_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(double *)(data1+i*stride1);
}

static void ebf_tfunc_int16_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}
static void ebf_tfunc_int16_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(int16_t *)(data1+i*stride1);
}


static void ebf_tfunc_int8_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}
static void ebf_tfunc_int8_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(int8_t *)(data1+i*stride1);
}


static void ebf_tfunc_uint8_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint8_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(uint8_t *)(data1+i*stride1);
}

static void ebf_tfunc_uint16_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint16_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(uint16_t *)(data1+i*stride1);
}


static void ebf_tfunc_uint32_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint32_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(uint32_t *)(data1+i*stride1);
}


static void ebf_tfunc_uint64_char(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(char *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_int32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int32_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_int64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int64_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_float32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(float *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_float64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(double *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_int16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int16_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_int8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(int8_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_uint8(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint8_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_uint16(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint16_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_uint32(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint32_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}
static void ebf_tfunc_uint64_uint64(char* data1,char* data2,size_t elements,size_t stride1,size_t stride2)
{
	size_t i;
	for(i=0;i<elements;++i)
		*(uint64_t *)(data2+i*stride2)=*(uint64_t *)(data1+i*stride1);
}



//------------------------------------------------------

int TypeS(const std::string& typestring)
{
	int etype = 8;
	if (typestring == "char")
	{
		etype = 1;
	}
	else if (typestring == "int32")
	{
		etype = 2;
	}
	else if (typestring == "int64")
	{
		etype = 3;
	}
	else if (typestring == "float32")
	{
		etype = 4;
	}
	else if (typestring == "float")
	{
		etype = 4;
	}
	else if (typestring == "float64")
	{
		etype = 5;
	}
	else if (typestring == "double")
	{
		etype = 5;
	}
	else if (typestring == "int16")
	{
		etype = 6;
	}
	else if (typestring == "int8")
	{
		etype = 9;
	}
	else if (typestring == "uint8")
	{
		etype = 10;
	}
	else if (typestring == "uint16")
	{
		etype = 11;
	}
	else if (typestring == "uint32")
	{
		etype = 12;
	}
	else if (typestring == "uint64")
	{
		etype = 13;
	}
//	else
//	{
//		throw std::logic_error("Ebf Error from Ebf_etype(): un recognized data type");
//	}
	return etype;

}





int64_t Ebf_GetElements(const std::string &filename1,const std::string &dataname1)
{
	int64_t ntot=0;
	int ecode=0;
	int64_t location=ebfc::EbfTable_Get(filename1.c_str(),dataname1.c_str(),&ecode);
	if (location>=0)
	{
		FILE* fp1=NULL;
		fp1=fopen(filename1.c_str(),"r");
		if(fp1!=NULL)
		{
			fseek(fp1,location,SEEK_SET);
			ebfc::EbfHeader ebfhm;
			ebfc::EbfHeader_read(&ebfhm,fp1);
			fclose(fp1);
			ntot=ebfc::EbfHeader_elements(&ebfhm);
		}
		else
		{
			std::cout<<filename1<<std::endl;
//			throw std::runtime_error("EBF Error from EbfFile::LoadMap: cannot find  above file");
		}
	}
	return ntot;
}





//-----------------------------------------------------------------------------
EbfFile::EbfFile()
{
	fp=NULL;
	debug=0;
}
EbfFile::EbfFile(const std::string& filename1, const std::string& dataname1, const std::string& mode1, int datatype, std::string dataunit, int rank, int64_t* dims)
{
	fp=NULL;
	debug=0;
	Open(filename1, dataname1, mode1, datatype, dataunit, rank, dims);
}
EbfFile::EbfFile(const EbfFile& other)
{
//		if(fp!=NULL)
//			throw std::range_error("Copying an inuse EbfFile not allowed");
	fp=NULL;
	debug=0;
}
EbfFile::~EbfFile()
{
	if(fp!=NULL)
	{
		Close();
		if ((mode==2)||(mode==3))
		{
			std::cout<<filename<<" "<<ebfh.headerpos<<" "<<ebfh.name<<" "<<ftell(fp)<<std::endl;
			std::cout<<"EBF WARNING: Efile Destructor is being used to close a file opened in write mode"<<std::endl;
			std::cout<<"It is preferable to use Efile.close() after a write operation"<<std::endl;
			throw std::logic_error("Efile Destructor is being used to close a file opened in write mode");
		}
	}
	else
	{
		Close();
	}
}
int64_t EbfFile::Index(int64_t i1,int64_t i2,int64_t i3)
{
	if(ebfh.rank==3)
		return (i1*ebfh.dim[1]+i2)*ebfh.dim[2]+i3;
	else
	{
		throw std::range_error("Ebf Error from index()::out of rank dim access");
	}
	return -1;
}
int64_t EbfFile::Index(int64_t i1,int64_t i2)
{
	if(ebfh.rank==2)
		return (i1*ebfh.dim[1]+i2);
	else
	{
		throw std::range_error("Ebf Error from index()::out of rank dim access");
	}
	return -1;
}



void EbfFile::Open(const std::string& filename1, const std::string& dataname1,
		std::string mode1, int datatype,std::string dataunit, int rank, int64_t *dim1)
{
	std::string sdef="";
	if (debug)
		std::cout << "Efile::open() " << filename1 + ":" + dataname1 << " mode="
				<< mode1 << std::endl;
	ecode=0;
	Close();
	filename = filename1;
	std::string mode2 = mode1;
	if (mode1 == "r")
	{
		mode = 1;
		mode2="rb";
		givenmode=1;
	}
	else if (mode1 == "w")
	{
		ecode=ebfc::EbfTable_Init(filename.c_str());
		if(ecode!=0)
			throw std::runtime_error("EBF Error from ebfc::EbfTable_Init()");
		mode = 3;
		mode2="rb+";
		givenmode=2;
	}
	else if (mode1 == "a")
	{
		mode = 3;
		mode2 = "rb+";
		givenmode=3;
	}
	else
	{
		throw std::logic_error(
				"EBF Error from Efile::open(): Mode must be r w or a here it is "
						+ mode1);
	}

	int64_t location=-1;
	if ((mode == 1) || (mode == 3))
		location = ebfc::EbfTable_Get(filename.c_str(), dataname1.c_str(),&ecode);

	if(ecode==0)
		fp = fopen(filename.c_str(), mode2.c_str());
	else
	{
		throw std::logic_error(
				"EBF Error from EBF::open(): probably not an ebf file");
	}

	if ((fp != NULL)&&(ecode==0))
	{
		if (mode == 3)
		{
			fseek(fp, 0, SEEK_END);
			std::vector<int64_t> dims1;
			if((dim1 == NULL)||(rank ==0))
			{
				dims1.push_back(0);
			}
			else
			{
				for(int i=0;i<rank;++i)
					dims1.push_back(dim1[i]);
			}

			ebfc::EbfHeader_set(&ebfh,dataname1.c_str(), datatype,
					ebfutils::TypeSize(datatype),dataunit.c_str(),sdef.c_str(),dims1.size(),&dims1[0]);
			if(ebfh.ecode!=0)
			{
				ecode=ebfh.ecode;
				printf("%s \n","EBF Error from ebfc::EbfHeader_set()");
			}
			if (location >= 0)
			{
				ecode=1;
				std::cout<<"EBF Error from EBF::open(): overwrite prevented as data object already exists OR not an ebf file "<<dataname1<<std::endl;
			}

			if(ecode==0)
				ebfc::EbfHeader_write(&ebfh,fp);

			if(ebfh.ecode!=0)
			{
				ecode=ebfh.ecode;
				printf("%s \n","EBF Error from ebfc::EbfHeader_set()");
			}
		}
		else if (mode == 1)
		{
			if (location >= 0)
			{
				fseek(fp, location, SEEK_SET);
				ebfc::EbfHeader_read(&ebfh,fp);
				if(ebfh.ecode!=0)
				{
					ecode=ebfh.ecode;
					printf("%s \n","EBF Error from ebfc::EbfHeader_read()");
				}
			}
			else
			{
				ecode=1;
				std::cout<<"EBF Error from EBF::open(): Key not found:"<<filename<<":"<<dataname1<<std::endl;
			}
		}
	}
	else
	{
		ecode=2;
		std::cout<<"EBF Error from EBF::open(): can't open file:"<<filename<<std::endl;
	}

	if(ecode!=0)
	{
		if(fp!=NULL)
		{
			fclose(fp);
			fp=NULL;
		}
		throw std::runtime_error("Ebf Error from Efile.open()");
	}
	if (debug)
		std::cout << "Efile::open() done" << std::endl;

}


void EbfFile::SaveTo(const std::string &filename1,const std::string& mode1)
{
	int ecode1=1;
	if (givenmode == 2)
	{
		Close();
		ecode1=ebfc::Ebf_Copy(filename.c_str(),filename1.c_str(),mode1.c_str(),"");
		if((ebfutils::FileExists(filename) == 1)&&(ecode1==0))
			remove(filename.c_str());
		if(ecode1!=0)
		{
			ecode=ecode1;
			throw std::runtime_error("EBF Error from EBF::saveTo(): while doing Ebf_Copy() ");
		}

	}
	else
	{
		Close();
		ecode=ecode1;
		throw std::runtime_error("EBF Error from EBF::saveTo():file must be opened in w mode for this operation");
	}

}


void EbfFile::Remove()
{
	if (givenmode == 2)
	{
		Close();
		if(ebfutils::FileExists(filename) == 1)
			remove(filename.c_str());
	}
	else
	{
		Close();
		ecode=1;
		throw std::runtime_error("EBF Error from EBF::saveTo():file must be opened in w mode for this operation");
	}

}


void EbfFile::Close()
{
	int myecode=0;
	int64_t mypos;
	if (fp != NULL)
	{
		if ((mode == 2) || (mode == 3))
		{
			if(ftell(fp)>ebfh.headerpos)
			{

			int64_t temp = 1;
			for (int i = 1; i < ebfh.rank; ++i)
				temp = temp * ebfh.dim[i];

			   mypos=ftell(fp);
		       if (temp*ebfh.datasize*ebfh.dim[0] != (mypos-ebfh.datapos))
		       {
		    	   if ((mypos-ebfh.datapos) == 0)
		    	   {
		    		   ebfh.rank=1;
		    		   ebfh.dim[0]=0;
		    	   }
		    	   else
		    	   {
		    		   if (((mypos-ebfh.datapos)%(temp*ebfh.datasize)) == 0)
		    		   {
		    			   ebfh.dim[0]=(mypos-ebfh.datapos)/(temp*ebfh.datasize);
		    		   }
		    		   else
		    		   {
		    			   int bufsize=(temp*ebfh.datasize)-(mypos-ebfh.datapos)%(temp*ebfh.datasize);
		    			   for(int i=0;i<bufsize;++i)
		    				   fwrite(&temp,1,1,fp);
							mypos=ftell(fp);
							ebfh.rank=1;
							ebfh.dim[0]=(mypos-ebfh.datapos)/ebfh.datasize;
							myecode=1;
		    		   }
		    	   }
					ebfh.capacity_=EbfHeader_elements(&ebfh) * ebfh.datasize;
					fseek(fp, ebfh.headerpos, SEEK_SET);
					ebfc::EbfHeader_write(&ebfh,fp);
		       }


			if ((mypos - ebfh.datapos)
					!= ebfh.capacity_)
				throw std::logic_error("Ebf error:efile_close() capacity_< total data size");


			}
		}

		if (fclose(fp) == 0)
			fp = NULL;
		else
		{
			  throw std::runtime_error("Ebf Error from Efile.close()-error closing file fclose()");
//			exit(1);
		}

		if (mode == 3)
		{
			ecode=ebfc::EbfTable_Put(filename.c_str(), ebfh.name,ebfh.headerpos);
		}
		if((ecode!=0)||(ebfh.ecode!=0))
			throw std::logic_error("Ebf Error from Efile.close()-either from Ebf-FileHT_put or ebfc::EbfHeader_write()");

		if(myecode!=0)
		{
			std::cout<<"Ebf error: Incorrect dimensions given, changing to 1 d and adjusting size"<<std::endl;
			throw std::logic_error("Ebf Error from Efile.close()");
		}


	}
	mode = 0;
	givenmode = 0;
}


void EbfFile::Seek(int64_t offset)
{
	if(mode==1)
	{
		if(offset>=0)
		{
			if(((offset)* ebfh.datasize+ebfh.datapos)>ebfh.datapos_end)
			{
				throw std::range_error("Ebf Error from Ebf::read()- trying to read past end of record");
			}
			else
			{
				fseek(fp, ebfh.datapos + offset * ebfh.datasize, SEEK_SET);
			}
		}
	}
	else
	{
		throw std::logic_error("Ebf Error from Ebf::position()- rewind allowed only in read mode");
	}
}

void EbfFile::ReadVoid(void* data,size_t data_size,size_t offset1)
{
	if (mode != 1)
	{
		throw std::logic_error("Ebf Error from Ebf::readBytes-Reading not on or status is zero");
	}

	if (offset1 >= 0)
		fseek(fp, ebfh.datapos + offset1, SEEK_SET);


	if (fread(data, ebfh.datasize * data_size, 1, fp)!= 1)
	{
		throw std::runtime_error("Ebf Error from Ebf::read()- fread");
	}
	if (ebfh.flag_swap == 1)
		ebfutils::SwapEndian(data,ebfh.datasize,data_size);
}



template<>
void EbfFile::ReadAsStringV(std::vector<std::string> &value1)
{

	std::stringstream sout;
	sout.setf(std::ios::scientific,std::ios::floatfield);
	sout.precision(7);

	switch(datatype())
	{
	case 1:
		{
			std::vector<char> buf1(value1.size());
			Read(&buf1[0], buf1.size());
			for (int64_t i = 0; i < int(value1.size()); ++i)
				sout <<buf1[i]<< " ";
		}
		break;
	case 2:
	{
		std::vector<int32_t> buf1(value1.size());
		Read(&buf1[0], buf1.size());
		for (int64_t i = 0; i < int(value1.size()); ++i)
			sout <<buf1[i]<< " ";
	}
		break;
	case 3:
	{
		std::vector<int64_t> buf1(value1.size());
		Read(&buf1[0], buf1.size());
		for (int64_t i = 0; i < int(value1.size()); ++i)
			sout <<buf1[i]<< " ";
	}
		break;
	case 4:
	{
		std::vector<efloat32> buf1(value1.size());
		Read(&buf1[0], buf1.size());
		for (int64_t i = 0; i < int(value1.size()); ++i)
			sout <<buf1[i]<< " ";
	}
		break;
	case 5:
	{
		sout.precision(15);
		std::vector<efloat64> buf1(value1.size());
		Read(&buf1[0], buf1.size());
		for (int64_t i = 0; i < int(value1.size()); ++i)
			sout <<buf1[i]<< " ";
	}
		break;
	case 6:
	{
		std::vector<int16_t> buf1(value1.size());
		Read(&buf1[0], buf1.size());
		for (int64_t i = 0; i < int(value1.size()); ++i)
			sout <<buf1[i]<< " ";
	}
		break;
	case 9:
	{
		std::vector<int8_t> buf1(value1.size());
		Read(&buf1[0], buf1.size());
		for (int64_t i = 0; i < int(value1.size()); ++i)
			sout <<int(buf1[i])<< " ";
	}
		break;
	case 10:
	{
		std::vector<uint8_t> buf1(value1.size());
		Read(&buf1[0], buf1.size());
		for (int64_t i = 0; i < int(value1.size()); ++i)
			sout <<int(buf1[i])<< " ";
	}
		break;
	case 11:
	{
		std::vector<uint16_t> buf1(value1.size());
		Read(&buf1[0], buf1.size());
		for (int64_t i = 0; i < int(value1.size()); ++i)
			sout <<buf1[i]<< " ";
	}
		break;
	case 12:
	{
		std::vector<uint32_t> buf1(value1.size());
		Read(&buf1[0], buf1.size());
		for (int64_t i = 0; i < int(value1.size()); ++i)
			sout <<buf1[i]<< " ";
	}
		break;
	case 13:
	{
		std::vector<uint64_t> buf1(value1.size());
		Read(&buf1[0], buf1.size());
		for (int64_t i = 0; i < int(value1.size()); ++i)
			sout <<buf1[i]<< " ";
	}
		break;
	default:
		throw std::runtime_error("Ebf Error from readstring-unrecognized datatype");
//		exit(1);
	}
	ebfutils::StringSplit(sout.str().c_str(), " ", value1);

}


void EbfFile::ReadGeneral(char* value1,size_t data_size,int typecode)
{
	static void (*EBFTFUNCARRAY[14][14])(char*, char*,size_t,size_t,size_t)={
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL,
			ebf_tfunc_char_char,
			ebf_tfunc_char_int32,
			ebf_tfunc_char_int64,
			ebf_tfunc_char_float32,
			ebf_tfunc_char_float64,
			ebf_tfunc_char_int16,
			NULL,
			NULL,
			ebf_tfunc_char_int8,
			ebf_tfunc_char_uint8,
			ebf_tfunc_char_uint16,
			ebf_tfunc_char_uint32,
			ebf_tfunc_char_uint64},
			{NULL,
			ebf_tfunc_int32_char,
			ebf_tfunc_int32_int32,
			ebf_tfunc_int32_int64,
			ebf_tfunc_int32_float32,
			ebf_tfunc_int32_float64,
			ebf_tfunc_int32_int16,
			NULL,
			NULL,
			ebf_tfunc_int32_int8,
			ebf_tfunc_int32_uint8,
			ebf_tfunc_int32_uint16,
			ebf_tfunc_int32_uint32,
			ebf_tfunc_int32_uint64},
			{NULL,
			ebf_tfunc_int64_char,
			ebf_tfunc_int64_int32,
			ebf_tfunc_int64_int64,
			ebf_tfunc_int64_float32,
			ebf_tfunc_int64_float64,
			ebf_tfunc_int64_int16,
			NULL,
			NULL,
			ebf_tfunc_int64_int8,
			ebf_tfunc_int64_uint8,
			ebf_tfunc_int64_uint16,
			ebf_tfunc_int64_uint32,
			ebf_tfunc_int64_uint64},
			{NULL,
			ebf_tfunc_float32_char,
			ebf_tfunc_float32_int32,
			ebf_tfunc_float32_int64,
			ebf_tfunc_float32_float32,
			ebf_tfunc_float32_float64,
			ebf_tfunc_float32_int16,
			NULL,
			NULL,
			ebf_tfunc_float32_int8,
			ebf_tfunc_float32_uint8,
			ebf_tfunc_float32_uint16,
			ebf_tfunc_float32_uint32,
			ebf_tfunc_float32_uint64},
			{NULL,
			ebf_tfunc_float64_char,
			ebf_tfunc_float64_int32,
			ebf_tfunc_float64_int64,
			ebf_tfunc_float64_float32,
			ebf_tfunc_float64_float64,
			ebf_tfunc_float64_int16,
			NULL,
			NULL,
			ebf_tfunc_float64_int8,
			ebf_tfunc_float64_uint8,
			ebf_tfunc_float64_uint16,
			ebf_tfunc_float64_uint32,
			ebf_tfunc_float64_uint64},
			{NULL,
			ebf_tfunc_int16_char,
			ebf_tfunc_int16_int32,
			ebf_tfunc_int16_int64,
			ebf_tfunc_int16_float32,
			ebf_tfunc_int16_float64,
			ebf_tfunc_int16_int16,
			NULL,
			NULL,
			ebf_tfunc_int16_int8,
			ebf_tfunc_int16_uint8,
			ebf_tfunc_int16_uint16,
			ebf_tfunc_int16_uint32,
			ebf_tfunc_int16_uint64},
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL,
			ebf_tfunc_int8_char,
			ebf_tfunc_int8_int32,
			ebf_tfunc_int8_int64,
			ebf_tfunc_int8_float32,
			ebf_tfunc_int8_float64,
			ebf_tfunc_int8_int16,
			NULL,
			NULL,
			ebf_tfunc_int8_int8,
			ebf_tfunc_int8_uint8,
			ebf_tfunc_int8_uint16,
			ebf_tfunc_int8_uint32,
			ebf_tfunc_int8_uint64},
			{NULL,
			ebf_tfunc_uint8_char,
			ebf_tfunc_uint8_int32,
			ebf_tfunc_uint8_int64,
			ebf_tfunc_uint8_float32,
			ebf_tfunc_uint8_float64,
			ebf_tfunc_uint8_int16,
			NULL,
			NULL,
			ebf_tfunc_uint8_int8,
			ebf_tfunc_uint8_uint8,
			ebf_tfunc_uint8_uint16,
			ebf_tfunc_uint8_uint32,
			ebf_tfunc_uint8_uint64},
			{NULL,
			ebf_tfunc_uint16_char,
			ebf_tfunc_uint16_int32,
			ebf_tfunc_uint16_int64,
			ebf_tfunc_uint16_float32,
			ebf_tfunc_uint16_float64,
			ebf_tfunc_uint16_int16,
			NULL,
			NULL,
			ebf_tfunc_uint16_int8,
			ebf_tfunc_uint16_uint8,
			ebf_tfunc_uint16_uint16,
			ebf_tfunc_uint16_uint32,
			ebf_tfunc_uint16_uint64},
			{NULL,
			ebf_tfunc_uint32_char,
			ebf_tfunc_uint32_int32,
			ebf_tfunc_uint32_int64,
			ebf_tfunc_uint32_float32,
			ebf_tfunc_uint32_float64,
			ebf_tfunc_uint32_int16,
			NULL,
			NULL,
			ebf_tfunc_uint32_int8,
			ebf_tfunc_uint32_uint8,
			ebf_tfunc_uint32_uint16,
			ebf_tfunc_uint32_uint32,
			ebf_tfunc_uint32_uint64},
			{NULL,
			ebf_tfunc_uint64_char,
			ebf_tfunc_uint64_int32,
			ebf_tfunc_uint64_int64,
			ebf_tfunc_uint64_float32,
			ebf_tfunc_uint64_float64,
			ebf_tfunc_uint64_int16,
			NULL,
			NULL,
			ebf_tfunc_uint64_int8,
			ebf_tfunc_uint64_uint8,
			ebf_tfunc_uint64_uint16,
			ebf_tfunc_uint64_uint32,
			ebf_tfunc_uint64_uint64}
	};

	if (mode != 1)
	{
		throw std::logic_error("Ebf Error from Ebf::read()-EBF error:Reading not on or status is zero");
	}


	if((int64_t(data_size)* ebfh.datasize+ftell(fp))>ebfh.datapos_end)
	{
		throw std::range_error("Ebf Error from Ebf::read()- trying to read past end of record");
	}

	if (ebfh.datatype == typecode)
	{
		if (fread(&value1[0], ebfh.datasize * data_size, 1, fp)!= 1)
		{
			throw std::runtime_error("Ebf Error from Ebf::read()- fread");
		}
		if (ebfh.flag_swap == 1)
			ebfutils::SwapEndian(&value1[0],ebfh.datasize,data_size);
	}
	else
	{
		int64_t ic = 0;
		int64_t size1 = sizeof(tbuf1.bufchar)/8;
		int64_t size = data_size;
		while (ic < size)
		{
			if ((ic + size1) > size)
				size1 = size - ic;
			if (fread(&tbuf1.bufchar[0], ebfh.datasize * size1, 1, fp)!= 1)
			{
				throw std::runtime_error("Ebf Error from Ebf::read()- fread");
			}
			if (ebfh.flag_swap == 1)
				ebfutils::SwapEndian(&tbuf1.bufchar[0],ebfh.datasize,size1);
			EBFTFUNCARRAY[ebfh.datatype][typecode](tbuf1.bufchar,value1+ic*ebfutils::TypeSize(typecode),size1,ebfutils::TypeSize(ebfh.datatype),ebfutils::TypeSize(typecode));
			ic = ic + size1;
		}
	}
}

void EbfFile::WriteGeneral(char* value1,size_t data_size,int typecode)
	{
	static void (*EBFTFUNCARRAY[14][14])(char*, char*,size_t,size_t,size_t)={
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL,
			ebf_tfunc_char_char,
			ebf_tfunc_char_int32,
			ebf_tfunc_char_int64,
			ebf_tfunc_char_float32,
			ebf_tfunc_char_float64,
			ebf_tfunc_char_int16,
			NULL,
			NULL,
			ebf_tfunc_char_int8,
			ebf_tfunc_char_uint8,
			ebf_tfunc_char_uint16,
			ebf_tfunc_char_uint32,
			ebf_tfunc_char_uint64},
			{NULL,
			ebf_tfunc_int32_char,
			ebf_tfunc_int32_int32,
			ebf_tfunc_int32_int64,
			ebf_tfunc_int32_float32,
			ebf_tfunc_int32_float64,
			ebf_tfunc_int32_int16,
			NULL,
			NULL,
			ebf_tfunc_int32_int8,
			ebf_tfunc_int32_uint8,
			ebf_tfunc_int32_uint16,
			ebf_tfunc_int32_uint32,
			ebf_tfunc_int32_uint64},
			{NULL,
			ebf_tfunc_int64_char,
			ebf_tfunc_int64_int32,
			ebf_tfunc_int64_int64,
			ebf_tfunc_int64_float32,
			ebf_tfunc_int64_float64,
			ebf_tfunc_int64_int16,
			NULL,
			NULL,
			ebf_tfunc_int64_int8,
			ebf_tfunc_int64_uint8,
			ebf_tfunc_int64_uint16,
			ebf_tfunc_int64_uint32,
			ebf_tfunc_int64_uint64},
			{NULL,
			ebf_tfunc_float32_char,
			ebf_tfunc_float32_int32,
			ebf_tfunc_float32_int64,
			ebf_tfunc_float32_float32,
			ebf_tfunc_float32_float64,
			ebf_tfunc_float32_int16,
			NULL,
			NULL,
			ebf_tfunc_float32_int8,
			ebf_tfunc_float32_uint8,
			ebf_tfunc_float32_uint16,
			ebf_tfunc_float32_uint32,
			ebf_tfunc_float32_uint64},
			{NULL,
			ebf_tfunc_float64_char,
			ebf_tfunc_float64_int32,
			ebf_tfunc_float64_int64,
			ebf_tfunc_float64_float32,
			ebf_tfunc_float64_float64,
			ebf_tfunc_float64_int16,
			NULL,
			NULL,
			ebf_tfunc_float64_int8,
			ebf_tfunc_float64_uint8,
			ebf_tfunc_float64_uint16,
			ebf_tfunc_float64_uint32,
			ebf_tfunc_float64_uint64},
			{NULL,
			ebf_tfunc_int16_char,
			ebf_tfunc_int16_int32,
			ebf_tfunc_int16_int64,
			ebf_tfunc_int16_float32,
			ebf_tfunc_int16_float64,
			ebf_tfunc_int16_int16,
			NULL,
			NULL,
			ebf_tfunc_int16_int8,
			ebf_tfunc_int16_uint8,
			ebf_tfunc_int16_uint16,
			ebf_tfunc_int16_uint32,
			ebf_tfunc_int16_uint64},
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL, NULL, NULL, NULL, NULL,},
			{NULL,
			ebf_tfunc_int8_char,
			ebf_tfunc_int8_int32,
			ebf_tfunc_int8_int64,
			ebf_tfunc_int8_float32,
			ebf_tfunc_int8_float64,
			ebf_tfunc_int8_int16,
			NULL,
			NULL,
			ebf_tfunc_int8_int8,
			ebf_tfunc_int8_uint8,
			ebf_tfunc_int8_uint16,
			ebf_tfunc_int8_uint32,
			ebf_tfunc_int8_uint64},
			{NULL,
			ebf_tfunc_uint8_char,
			ebf_tfunc_uint8_int32,
			ebf_tfunc_uint8_int64,
			ebf_tfunc_uint8_float32,
			ebf_tfunc_uint8_float64,
			ebf_tfunc_uint8_int16,
			NULL,
			NULL,
			ebf_tfunc_uint8_int8,
			ebf_tfunc_uint8_uint8,
			ebf_tfunc_uint8_uint16,
			ebf_tfunc_uint8_uint32,
			ebf_tfunc_uint8_uint64},
			{NULL,
			ebf_tfunc_uint16_char,
			ebf_tfunc_uint16_int32,
			ebf_tfunc_uint16_int64,
			ebf_tfunc_uint16_float32,
			ebf_tfunc_uint16_float64,
			ebf_tfunc_uint16_int16,
			NULL,
			NULL,
			ebf_tfunc_uint16_int8,
			ebf_tfunc_uint16_uint8,
			ebf_tfunc_uint16_uint16,
			ebf_tfunc_uint16_uint32,
			ebf_tfunc_uint16_uint64},
			{NULL,
			ebf_tfunc_uint32_char,
			ebf_tfunc_uint32_int32,
			ebf_tfunc_uint32_int64,
			ebf_tfunc_uint32_float32,
			ebf_tfunc_uint32_float64,
			ebf_tfunc_uint32_int16,
			NULL,
			NULL,
			ebf_tfunc_uint32_int8,
			ebf_tfunc_uint32_uint8,
			ebf_tfunc_uint32_uint16,
			ebf_tfunc_uint32_uint32,
			ebf_tfunc_uint32_uint64},
			{NULL,
			ebf_tfunc_uint64_char,
			ebf_tfunc_uint64_int32,
			ebf_tfunc_uint64_int64,
			ebf_tfunc_uint64_float32,
			ebf_tfunc_uint64_float64,
			ebf_tfunc_uint64_int16,
			NULL,
			NULL,
			ebf_tfunc_uint64_int8,
			ebf_tfunc_uint64_uint8,
			ebf_tfunc_uint64_uint16,
			ebf_tfunc_uint64_uint32,
			ebf_tfunc_uint64_uint64}
	};

		if((mode!=2)&&(mode!=3))
		{
			throw std::logic_error("Ebf Error from Ebf::write()-EBF error: Writing not on or status is zero write blocked");
		}

		if (ebfh.datatype == typecode)
		{
			if(fwrite(&value1[0],ebfh.datasize * data_size,1,fp)!=1)
			{
				throw std::runtime_error("Ebf Error from Ebf::write()- fwrite error");
			}
		}
		else
		{
			int64_t ic = 0;
			int64_t size1 = sizeof(tbuf1.bufchar)/8;
			int64_t size = data_size;
			while (ic < size)
			{
				if ((ic + size1) > size)
					size1 = size - ic;
				EBFTFUNCARRAY[typecode][ebfh.datatype](value1+ic*ebfutils::TypeSize(typecode),tbuf1.bufchar,size1,ebfutils::TypeSize(typecode),ebfutils::TypeSize(ebfh.datatype));
				if (fwrite(&tbuf1.bufchar[0], ebfh.datasize * size1, 1, fp)!= 1)
				{
					throw std::logic_error("Ebf Error from Ebf::write()-status not 1, file may be corrupted");
				}
				ic = ic + size1;
			}
		}
	}



	void EbfFile::Write(char* value1, int64_t data_size)
	{
		WriteGeneral((char *)value1,data_size,1);
	}
	void EbfFile::Write(int32_t* value1, int64_t data_size)
	{
		WriteGeneral((char *)value1,data_size,2);
	}
	void EbfFile::Write(int64_t* value1, int64_t data_size)
	{
		WriteGeneral((char *)value1,data_size,3);
	}
	void EbfFile::Write(efloat32* value1, int64_t data_size)
	{
		WriteGeneral((char *)value1,data_size,4);
	}
	void EbfFile::Write(efloat64* value1, int64_t data_size)
	{
		WriteGeneral((char *)value1,data_size,5);
	}
	void EbfFile::Write(int16_t* value1, int64_t data_size)
	{
		WriteGeneral((char *)value1,data_size,6);
	}
	void EbfFile::Write(int8_t* value1, int64_t data_size)
	{
		WriteGeneral((char *)value1,data_size,9);
	}
	void EbfFile::Write(uint8_t* value1, int64_t data_size)
	{
		WriteGeneral((char *)value1,data_size,10);
	}
	void EbfFile::Write(uint16_t* value1, int64_t data_size)
	{
		WriteGeneral((char *)value1,data_size,11);
	}
	void EbfFile::Write(uint32_t* value1, int64_t data_size)
	{
		WriteGeneral((char *)value1,data_size,12);
	}
	void EbfFile::Write(uint64_t* value1, int64_t data_size)
	{
		WriteGeneral((char *)value1,data_size,13);
	}


	void EbfFile::Read(std::string *value1, size_t data_size = 1)
	{
		throw std::logic_error("Ebf Error from readstring():-Forbidden to use string data type");
	}
	void EbfFile::Read(char *value1, size_t data_size = 1)
	{
		ReadGeneral((char *)value1,data_size,1);
	}
	void EbfFile::Read(int32_t *value1, size_t data_size = 1)
	{
		ReadGeneral((char *)value1,data_size,2);
	}
	void EbfFile::Read(int64_t *value1, size_t data_size = 1)
	{
		ReadGeneral((char *)value1,data_size,3);
	}
	void EbfFile::Read(efloat32 *value1, size_t data_size = 1)
	{
		ReadGeneral((char *)value1,data_size,4);
	}
	void EbfFile::Read(efloat64 *value1, size_t data_size = 1)
	{
		ReadGeneral((char *)value1,data_size,5);
	}
	void EbfFile::Read(int16_t *value1, size_t data_size = 1)
	{
		ReadGeneral((char *)value1,data_size,6);
	}
	void EbfFile::Read(int8_t *value1, size_t data_size = 1)
	{
		ReadGeneral((char *)value1,data_size,9);
	}
	void EbfFile::Read(uint8_t *value1, size_t data_size = 1)
	{
		ReadGeneral((char *)value1,data_size,10);
	}
	void EbfFile::Read(uint16_t *value1, size_t data_size = 1)
	{
		ReadGeneral((char *)value1,data_size,11);
	}
	void EbfFile::Read(uint32_t *value1, size_t data_size = 1)
	{
		ReadGeneral((char *)value1,data_size,12);
	}
	void EbfFile::Read(uint64_t *value1, size_t data_size = 1)
	{
		ReadGeneral((char *)value1,data_size,13);
	}



int ReadString(const std::string& filename, const std::string& dataname,std::string &mystr)
{
	std::vector<char> x;
	EbfFile ebf1(filename, dataname);
	x.resize(ebf1.elements()+1);
	if(ebf1.ecode==0)
		ebf1.Read(&x[0],ebf1.elements());
	ebf1.Close();
	x[x.size()-1]=0;
	mystr=(&x[0]);
	return ebf1.ecode;
}

int ReadString(const std::string& filename, const std::string& dataname,std::vector<std::string> &mystr)
{
	std::vector<char> x;
	EbfFile ebf1(filename, dataname);
	x.resize(ebf1.elements());
	if(ebf1.ecode==0)
		ebf1.Read(&x[0],ebf1.elements());

	std::string temp;
	mystr.clear();
	for(int64_t i=0;i<ebf1.dim(1);++i)
	{
		temp.assign(&x[i*ebf1.dim(1)],std::min(int64_t(strlen(&x[i*ebf1.dim(1)])),ebf1.dim(1)));
		mystr.push_back(temp);
	}
	ebf1.Close();

	return ebf1.ecode;
}


int WriteString(const std::string& filename, const std::string& dataname, const std::string &data, const std::string& mode)
{
	std::vector<char> x(data.size());
	data.copy(&x[0],data.size());
	if(x.size() >0)
		return Write(filename,dataname,&x[0],mode,"",x.size());
	else
		return 1;
}

int WriteString(const std::string& filename, const std::string& dataname, const std::vector<std::string> &data, const std::string& mode)
{
	size_t maxlen=0;
	for(size_t i=0;i<data.size();++i)
	{
		if(maxlen<data[i].size())
			maxlen=data[i].size();
	}
	std::vector<char> x(data.size()*maxlen,0);
	for(size_t i=0;i<data.size();++i)
	{
		data[i].copy(&x[i*maxlen],data[i].size());
	}

	if(x.size() >0)
		return Write(filename,dataname,&x[0],mode,"",x.size()/maxlen,maxlen);
	else
		return 1;

}



int ContainsKey(const std::string& infile, const std::string& dataname)
{
	ebfc::EbfDataInfo dinfo;
	return ebfc::Ebf_ContainsKey_Info(infile.c_str(),dataname.c_str(),&dinfo);
}
int ContainsKey(const std::string& infile, const std::string& dataname,ebfc::EbfDataInfo &dinfo)
{
	return ebfc::Ebf_ContainsKey_Info(infile.c_str(),dataname.c_str(),&dinfo);
}


int Rename(const std::string  &filename,const  std::string &oldkey,const std::string  &newkey)
{
	return  ebfc::Ebf_Rename(filename.c_str(),oldkey.c_str(),newkey.c_str());
}

//---------------------------------------------------------------



}
