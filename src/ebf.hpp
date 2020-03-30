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


#ifndef EBF_H_
#define EBF_H_
#include <vector>
#include <string>
#include "ebftable.hpp"


namespace ebf
{




namespace ebfutils
{
void Ebf_Error(const std::string& mystr);
bool FileExists(const std::string& filename);
std::string EbfCodeToName(int x);
void StringSplit(const std::string &s,const char* delimiters,std::vector<std::string> &sv);
void SwapEndian(void* addr, int datasize,int64_t nmemb=1);
int TypeSize(int x);
}


template<class T>
int TypeT(void) {return -1;}
template<> int TypeT<char>(void);
template<> int TypeT<int32_t> (void);
template<> int TypeT<int64_t> (void);
template<> int TypeT<efloat32> (void);
template<> int TypeT<efloat64> (void);
template<> int TypeT<int16_t> (void);
template<> int TypeT<std::string> (void);
template<> int TypeT<int8_t> (void);
template<> int TypeT<uint8_t> (void);
template<> int TypeT<uint16_t> (void);
template<> int TypeT<uint32_t> (void);
template<> int TypeT<uint64_t> (void);

template<class T>
int TypeV(const T &x) {return TypeT<T>();}



int TypeS(const std::string &typestring);


/**
 *
 * @param filename1  name of file  e.g "check.ebf"
 * @param dataname1  data tag name within the file e.g "/x"
 * @return number of GetElements, if return value is zero it indicates error
 * it could be non existent file, non existent item, or zero size data object
 */
// int64_t Ebf_GetElements(const std::string &filename1,const std::string &dataname1);

//inline bool Ebf_ContainsKey(const std::string &filename1,const std::string &dataname1)
//{
//	int ecode=0;
//	int64_t loc=ebfc::EbfTable_Get(filename1.c_str(),dataname1.c_str(),&ecode);
//	if(ecode!=0)
//		throw std::logic_error("Ebf Error: while running EbfTable_Get()");
//	return loc;
//}

// void Ebf_GetTagNames(const char* filename, std::vector<std::string> &tagnames,int *ecode);
void Info(const std::string& filename1);

typedef union
{
	char bufchar[8000];
	int8_t bufint8[8000];
	int16_t bufint16[4000];
	int32_t bufint32[2000];
	int64_t bufint64[1000];
	uint8_t bufuint8[8000];
	uint16_t bufuint16[4000];
	uint32_t bufuint32[2000];
	uint64_t bufuint64[1000];
	efloat32 buffloat32[2000];
	efloat64 buffloat64[1000];
} EbfFileBuf;

class EbfFile
{
private:
public:
	EbfFile();
	EbfFile(const std::string& filename1, const std::string& dataname1, const std::string& mode1="r", int datatype=0, std::string dataunit="", int rank=1, int64_t* dims=NULL);
	EbfFile(const EbfFile& other);
	~EbfFile();
	int debug,ecode;
private:
	FILE* fp;
	ebfc::EbfHeader ebfh;
	int mode,givenmode;
	std::string filename;
//	std::vector<char> tbuf;
	EbfFileBuf tbuf1;
//	std::vector<int64_t> patternData;
//	int64_t ncur;
public:
	// main reading/writing functions
	/**
	 * Open a file and initialize file buffers to perform Read Write operations
	 * also positions the file to the location of the data.
	 * If a file is already Open from previous operation it is closed before opening a new buffer.
	 * @param filename1
	 * @param dataname1
	 * @param mode1     "w", "a" or "r"
	 * @param datatype the desired destination datatype
	 * @param dataunit units as a string
	 * @param rank     the rank or dimensionality of data
	 * @param dims     the shape or dimensionality of arrays being written.
	 * Depending upon the number of elements written, the value of dims[0] is adjusted such
	 * that the product dim[0]..dim[rank-1]= no of elements written. This is done the ebf file is closed.
	 */
	void Open(const std::string& filename1, const std::string& dataname1, std::string mode1="r", int datatype=0, std::string dataunit="", int rank=1, int64_t* dims=NULL);
	void SaveTo(const std::string &filename1, const std::string& mode1);
	void Remove();
	void Close();
	// enquire header properties
	inline int64_t headerpos()	  { return ebfh.headerpos;}
	inline int64_t datapos()	  { return ebfh.datapos;}
	inline int64_t elements()	  { return EbfHeader_elements(&ebfh);}
	inline int datatype()	  {	return ebfh.datatype;}
	inline int datasize()	  {	return ebfh.datasize;}
	inline int capacity()	  {	return ebfh.capacity_;}
	inline int rank()	      {	return ebfh.rank;}
	inline  std::string getFileName() {	return filename;}
	inline  std::string getDataName() {	return ebfh.name;}
	inline int64_t dim(int i)	  {	return ebfh.dim[i];}
	inline std::string unit() { return ebfh.unit;}
	inline std::string sdef()  { return ebfh.sdef; }
    int64_t Index(int64_t i1,int64_t i2,int64_t i3);
    int64_t Index(int64_t i1,int64_t i2);
    // set Tell	get status
	void Seek(int64_t offset);
	int64_t Tell()	          { return (ftell(fp)-ebfh.datapos)/ebfh.datasize;}
	// Read and Write templates
	void Write(char* value1, int64_t data_size);
	void Write(int32_t* value1, int64_t data_size);
	void Write(int64_t* value1, int64_t data_size);
	void Write(efloat32* value1, int64_t data_size);
	void Write(efloat64* value1, int64_t data_size);
	void Write(int16_t* value1, int64_t data_size);
	void Write(int8_t* value1, int64_t data_size);
	void Write(uint8_t* value1, int64_t data_size);
	void Write(uint16_t* value1, int64_t data_size);
	void Write(uint32_t* value1, int64_t data_size);
	void Write(uint64_t* value1, int64_t data_size);
	void Read(std::string *value1, size_t data_size);
	void Read(char *value1, size_t data_size);
	void Read(int32_t *value1, size_t data_size);
	void Read(int64_t *value1, size_t data_size);
	void Read(efloat32 *value1, size_t data_size);
	void Read(efloat64 *value1, size_t data_size);
	void Read(int16_t *value1, size_t data_size);
	void Read(int8_t *value1, size_t data_size);
	void Read(uint8_t *value1, size_t data_size);
	void Read(uint16_t *value1, size_t data_size);
	void Read(uint32_t *value1, size_t data_size);
	void Read(uint64_t *value1, size_t data_size);
// Read raw bytes correctly byte aligned
// Note: size and offset are in units of bytes
	void ReadVoid(void* data,size_t data_size,size_t offset1=-1);
	template<class T1>
	void ReadAsStringV(std::vector<T1> &value1)
	{
		ebfutils::Ebf_Error("Ebf Error from readstring():-Forbidden to use non string data type");
	}
private:
	void ReadGeneral(char* value1,size_t data_size,int typecode);
	void WriteGeneral(char* value1,size_t data_size,int typecode);
//	size_t AutoSwappedFread(void *ptr, int64_t *pattern, size_t nmemb, int flag_swap,
//			FILE *stream);

};

template<>
void EbfFile::ReadAsStringV(std::vector<std::string> &value1);



template<class T1, class T2> int WriteAs(const std::string& filename,
		const std::string& dataname, T2 * value1, const std::string& mode, std::string dataunit="",int64_t dim1=1,
		int64_t dim2 = 0, int64_t dim3 = 0, int64_t dim4 = 0)
{
	int64_t dims[8];
	int rank=1;
	dims[0]=dim1;
	if(dim2>0)	{		dims[rank]=dim2;		rank++;	}
	if(dim3>0)	{		dims[rank]=dim3;		rank++;	}
	if(dim4>0)	{		dims[rank]=dim4;		rank++;	}
//	std::cout<<"Writing: "<<filename+":"+dataname<<" tc="<<TypeT<T2>()<<std::endl;
	EbfFile ebf1(filename, dataname, mode, TypeT<T1>(), dataunit,rank,dims);
	ebf1.Write(&value1[0],ebf1.elements());
	ebf1.Close();
	return ebf1.ecode;
//	std::cout<<"Written: "<<filename+":"+dataname<<" tc="<<TypeT<T2>()<<std::endl;
}

template<class T2> int Write(const std::string& filename,
		const std::string& dataname, T2 * value1, const std::string& mode, std::string dataunit="", int64_t dim1=1,
		int64_t dim2 = 0, int64_t dim3 = 0, int64_t dim4 = 0)
{
	int64_t dims[8];
	int rank=1;
//	std::cout<<"Writing: "<<filename+":"+dataname<<" tc="<<TypeT<T2>()<<std::endl;
	dims[0]=dim1;
	if(dim2>0)	{		dims[rank]=dim2;		rank++;	}
	if(dim3>0)	{		dims[rank]=dim3;		rank++;	}
	if(dim4>0)	{		dims[rank]=dim4;		rank++;	}
	EbfFile ebf1(filename, dataname, mode, TypeT<T2>(),dataunit,rank,&dims[0]);
	if(ebf1.ecode==0)
		ebf1.Write(&value1[0],ebf1.elements());
	ebf1.Close();
	return ebf1.ecode;
//	std::cout<<"Written: "<<filename+":"+dataname<<" tc="<<TypeT<T2>()<<std::endl;
}

template<class T1> int Read(const std::string& filename, const std::string& dataname,	std::vector<T1> &value1)
{
//	std::cout<<"Reading: "<<filename+":"+dataname<<std::endl;
	EbfFile ebf1(filename, dataname);
	int64_t ntot=ebf1.elements();
	if(int64_t(value1.size())!=ntot)
		value1.resize(ntot);
	if(ebf1.ecode==0)
		ebf1.Read(&value1[0],ntot);
	ebf1.Close();
	return ebf1.ecode;
//	std::cout<<"Read: "<<filename+":"+dataname<<std::endl;
}


template<class T1> int Read(const std::string& filename, const std::string& dataname,T1* value1,int64_t ntot=1, int64_t offset1=-1)
{
//	std::cout<<"Reading: "<<filename+":"+dataname<<std::endl;
	EbfFile ebf1(filename, dataname);
	if(ntot>ebf1.elements())
	{
		ebf1.Close();
		ebfutils::Ebf_Error("Ebf Error from ebfread()::supplied size is larger than size of data on file");
	}
	if(offset1>=0)
		ebf1.Seek(offset1);
	if(ebf1.ecode==0)
		ebf1.Read(&value1[0],ntot);
	ebf1.Close();
	return ebf1.ecode;
//	std::cout<<"Read: "<<filename+":"+dataname<<std::endl;
}


int ReadString(const std::string& filename, const std::string& dataname,std::string &mystr);
int ReadString(const std::string& filename, const std::string& dataname,std::vector<std::string> &mystr);
int WriteString(const std::string& filename, const std::string& dataname, const std::string &data, const std::string& mode);
int WriteString(const std::string& filename, const std::string& dataname, const std::vector<std::string> &data, const std::string& mode);


//void Copy(const std::string& infile,const std::string& outfile,const std::string& mode1);
//void Copy(const std::string& infile,const std::string& outfile,const std::string& mode1,std::string dataname="");

/**
 * Copy data items from one file to another
 *
 *  @param infile
 *  @param outfile
 *  @param mode     "w" or "a"
 *  @param dataname the desired list of tagnames to transfer, blank string to trnasfer all items
 */
int Copy(const std::string& infile, const std::string& outfile,const std::string& mode="a", const std::string& dataname="");

int ContainsKey(const std::string& infile, const std::string& dataname);
int ContainsKey(const std::string& infile, const std::string& dataname,ebfc::EbfDataInfo &dinfo);
int Rename(const std::string  &filename,const  std::string &oldkey,const std::string  &newkey);

using ebfc::EbfDataInfo;

//--------------------------------------------------------------------------------

}

//namespace ebf
//{
//using ebf_impl1::EbfFile;
//using ebf_impl1::ebfc::EbfDataInfo;
//using ebf_impl1::ContainsKey;
//using ebf_impl1::Copy;
//using ebf_impl1::Info;
//using ebf_impl1::Read;
//using ebf_impl1::Write;
//using ebf_impl1::WriteAs;
//using ebf_impl1::ReadString;
//using ebf_impl1::WriteString;
//using ebf_impl1::TypeT;
//using ebf_impl1::TypeV;
//using ebf_impl1::TypeS;
//using ebf_impl1::Rename;
//}

#endif /* EBF_H_ */










//template<class T1> void ebfralloc(const std::string& filename, const std::string& dataname)
//{
//	T1 *value1=NULL;
//	try
//	{
//		EbfFile ebf1(filename, dataname);
//		int64_t ntot=ebf1.GetElements();
//		T1* value1= new T1[ntot];
//		ebf1.Read(&value1[0],ntot);
//	}
//	catch (std::bad_alloc& ba)
//	{
//	  std::cout << "Error allocating memory."<<ba.what() << std::endl;
//	  exit(1);
//	}
//	return value1;
//}
