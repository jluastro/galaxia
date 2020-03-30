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


#ifndef EBFHT_H_
#define EBFHT_H_


#define EBFLANGOPTCPP

/*--------------------	--------------------------------------------*/
/* using cstdint cinttypes gives problem only works for -std=c++0x */
/*
#include <stdint.h>
#include <inttypes.h>
*/


#ifdef EBFLANGOPTCPP
#include <stdint.h>
#include <inttypes.h>
#include<cstdio>
namespace ebf{
#else
#include <stdint.h>
#include <inttypes.h>
#include<stdio.h>
#ifdef __cplusplus
 extern "C" {
#endif
#endif

typedef float efloat32; /* 64 bit integer unsigned */
typedef double efloat64; /* 64 bit integer unsigned */

#ifdef EBFLANGOPTCPP
namespace ebfc{
#endif




#define EBF_KEYNAME_MAXSIZE 512
#define EBF_SDEF_MAXSIZE 2048
#define EBF_RANK_MAXSIZE 8
#define EBF_VERSION_MAJOR 0
#define EBF_VERSION_MINOR 0
#define EBF_VERSION_REVISION 0




/**
 * A structure to hold information about a data object and return it to user
 * ecode     - if non zero then either data not present or problem loading data
 * location  - position in file pointing to start of the data object i.e, headerpos
 * datapos   - position in file of the pointing to start of the raw data
 * datatype  - data type of the object
 * datasize  - size in bytes of each data item corresponding ot the datatype i.e., 4 for float
 * elements  - total number of elements or items in the data object
 * rank      - number of dimensions of the data
 * dim       - an array specifying the dimensions
*/
typedef struct
{
    int ecode,datatype,datasize,rank;
    int64_t headerpos,datapos,elements;
    char unit[EBF_KEYNAME_MAXSIZE];
    int64_t dim[EBF_RANK_MAXSIZE];
} EbfDataInfo;


/**
 * Usage internal: as structure used for reading ebf headers
 */
typedef struct
{
	char sig[8];
    char version[4];
    char flags[4];
    int  ecode;
    char name[EBF_KEYNAME_MAXSIZE];
    char unit[EBF_KEYNAME_MAXSIZE];
    char sdef[EBF_SDEF_MAXSIZE];
    int32_t endiantest,extrasize;
    int32_t headersize,datatype,datasize,rank;
    int flag_swap;
    int64_t capacity_;
    int64_t dim[EBF_RANK_MAXSIZE];
    int64_t headerpos,datapos,datapos_end;
} EbfHeader;

int64_t EbfHeader_elements(EbfHeader *ebfh );
void EbfHeader_set(EbfHeader* ebfh, const char* dataname1, int datatype1, int datasize1,const char* dataunit,const char* sdef, int rank, const int64_t* dims1);
void EbfHeader_read(EbfHeader *ebfh,FILE* fp);
void EbfHeader_write(EbfHeader *ebfh,FILE* fp);



int EbfTable_Init(const char*  filename);
int EbfTable_InitSwap(const char*  filename);
int EbfTable_Put(const char*  filename,const char*  key, int64_t value);

/**
 * A function to check existence of an item in a file and get its location
 * @param filename  name of file  e.g "check.ebf"
 * @param key  data object name within the file e.g "/x"
 * @param ecode  non zero if there was an error
 * @return location pointing to data object in file a value <0 means data object not found
 */
int64_t EbfTable_Get(const char*  filename,const char*  key, int *ecode);
int EbfTable_Remove(const char*  filename,const char*  key);
int Ebf_Rename(const char*  filename,const char*  oldkey,const char*  newkey);

/**
 *  A function to copy a data object from one ebf file to another
 * @param infile source file from which to copy data
 * @param outfile destination file to which data is to be copied
 * @param mode1 I/O mode "w" or "a" , fresh write or append respectively
 * @param dataname the name of data to be copied can be comma seprated values also
 * @return an error code to test if operation had any error. Non zero on error.
 *        If blank string then all data items in infile are copied.
 */
int Ebf_Copy(const char* infile, const char* outfile, const char* mode1,const char* dataname);

/*
 * A function to check existence of a data object in a file
 * @param filename1  name of file  e.g "check.ebf"
 * @param dataname1  data object name within the file e.g "/x"
 * @return true if item exists or else false
 * @param ecode an error code to test if operation had any error. Non zero on error.
 */
/* int Ebf_ContainsKey(const char* filename1, const char* dataname1,int* ecode); */

int64_t Ebf_GetCheckSum(const char* filename,int *ecode);


int Ebf_Version(int i);


/**
 * A function to check existence of an item in a file and get its location.
 * Also gets the EbfDataInfo structure, which can be queried for number of elements
 * data units and existence of data and so on.
 * @param filename1  name of file  e.g "check.ebf"
 * @param dataname1  data tag name within the file e.g "/x"
 * @param dinfo a structure in which the information about data will be put
 * @return an integer 1 if data item is present, if not 1 then it means that either data
 *         item is not present or file is not present or problem loading the data item.
 * Among other things dinfo can be queried for
 *         - dinfo.elements
 *         - dinfo.unit
 *         .
 */
int Ebf_ContainsKey_Info(const char * filename1,const char* dataname1,EbfDataInfo* dinfo);


void EbfTable_Print(const char* filename);

#ifdef EBFLANGOPTCPP
}
}
#else
#ifdef __cplusplus
}
#endif
#endif /* namespace end */


#endif /* EBFHT_H_ */


