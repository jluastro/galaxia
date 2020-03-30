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

#include "ebftable.hpp"

#ifdef EBFLANGOPTCPP
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <cstring>
#include<cstdlib>
#include<cctype>
namespace ebf{
namespace ebfc{
#else
#include <string.h>
#include<stdlib.h>
#include<ctype.h>
#endif

/*
#ifdef __cplusplus
#include <cstring>
#include<cstdlib>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include<cctype>
#include "ebftable.h"
namespace ebf_impl1{
namespace ebfc{
#else
#include <string.h>
#include<stdlib.h>
#include<ctype.h>
#include "ebftable.h"
#endif
*/

int Ebf_Version(int i)
{
	int version[3]={0,0,0};
	if((i<0)||(i>2))
	{
		return -1;
	}
	return version[i];
}

void ebfutils_SwapEndian(void* addr, int datasize, int64_t nmemb)
{
	int64_t pattern[3];
	int64_t i, j = 0, k;
	char c;
	pattern[0] = datasize;
	pattern[1] = nmemb;
	pattern[2] = 0;

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


/*--------------------------------------------------------------------*/

static char * ebfstrdup(const char *src)
{
	char *copy = NULL;
	if (src != NULL)
	{
		copy = (char *) malloc(strlen(src) + 1);
		if (copy != NULL)
			strcpy(copy, src);
	}
	return copy;
}

static void ebfstrtrim(char*  dest,size_t num)
{
	char *temp=strchr(dest,' ');
	if((temp!=NULL)&&((temp-dest)<(int)num))
		*temp=0;
}

static int ebfstrcpy(char*  dest,size_t num,const char* src)
{
	if(num>strlen(src))
	{
		strcpy(dest,src);
		return 0;
	}
	else
	{
		strncpy(dest,src,num);
		dest[num-1]=0;
		printf("%s \n","Ebf error: in ebfstrcpy() src string too large");
		return 1;
	}
}


/*--------------------------------------------------------------------*/




int64_t EbfHeader_elements(EbfHeader* ebfh)
{
	int i = 0;
	int64_t blklen = 1;
/*	if((ebfh->ecode==0)&&(ebfh->rank <= EBF_RANK_MAXSIZE)) */
	if(ebfh->rank <= EBF_RANK_MAXSIZE)
	{
		for (i = 0; i < ebfh->rank; ++i)
			blklen *= (int64_t) (ebfh->dim[i]);
		return blklen;
	}
	else
	{
		printf("%s \n","Ebf error: in EbfHeader_elements()");
		return 0;
	}
}


void EbfHeader_print(EbfHeader* ebfh)
{
	int i = 0;
	printf(
			"%s \n",
			"-----------------------------------------------------------------");
	printf("%12s %12s", "name=", ebfh->name);
	printf("%12s %12d %12s %12d %12s %12d %12s ", "datatype=", ebfh->datatype,
			"datasize=", ebfh->datasize, "rank=", ebfh->rank, "dim=");

	printf("%s", "(");
	for (i = 0; i < ebfh->rank - 1; ++i)
		printf("%d ", (int) ebfh->dim[i]);
	printf("%s \n", ")");
}

void EbfHeader_set(EbfHeader* ebfh, const char* dataname1, int datatype1,
		int datasize1,const char* dataunit,const char* sdef, int rank, const int64_t* dims1)
{
	int i = 0;
	size_t j = 0;
	int extrasize;
	ebfh->ecode=0;
	ebfh->sig[0] = -118;
	ebfh->sig[1] = 69;
	ebfh->sig[2] = 66;
	ebfh->sig[3] = 70;
	ebfh->sig[4] = -82;
	ebfh->sig[5] = 43;
	ebfh->sig[6] = -81;
	ebfh->sig[7] = 10;
	ebfh->flag_swap=0;


	if(ebfstrcpy(ebfh->name, EBF_KEYNAME_MAXSIZE,dataname1)!=0)
		ebfh->ecode=307;
	if(ebfstrcpy(ebfh->unit, EBF_KEYNAME_MAXSIZE,dataunit)!=0)
		ebfh->ecode=307;
	if(ebfstrcpy(ebfh->sdef, EBF_SDEF_MAXSIZE,sdef)!=0)
		ebfh->ecode=307;
	for(j=0;j<(strlen(ebfh->name));++j)
		ebfh->name[j]=tolower(ebfh->name[j]);


	ebfh->endiantest = 1684234849;
	ebfh->datatype = datatype1;
	ebfh->datasize = datasize1;
	ebfh->rank = rank;
	ebfh->headerpos = 0;
	if(ebfh->rank > EBF_RANK_MAXSIZE)
	{
		ebfh->ecode=311;
		ebfh->rank=EBF_RANK_MAXSIZE;
	}

	for (i = 0; i < ebfh->rank; ++i)
		ebfh->dim[i] = dims1[i];

	ebfh->version[0] = 1;
	ebfh->version[1] = 1;
	ebfh->version[2] = 0;
	ebfh->version[3] = 0;
	ebfh->flags[0] = 0;
	ebfh->flags[1] = 0;
	ebfh->flags[2] = 0;
	ebfh->flags[3] = 0;

	if ((ebfh->rank == 1) && (ebfh->dim[0] == 1))
		extrasize = 16;
	else
		extrasize = 64;

	ebfh->capacity_=EbfHeader_elements(ebfh)*ebfh->datasize;
	ebfh->headersize = 56 + strlen(ebfh->name) + strlen(ebfh->unit)+ strlen(ebfh->sdef)
				+ extrasize + 8 * ebfh->rank;

}

static void EbfHeader_writeVer11(EbfHeader* ebfh, FILE* fp)
{
	int32_t unused4[10];
	int64_t unused8[10];
	int nwrit = 1;
	int32_t namesize = strlen(ebfh->name);
	int32_t unitsize = strlen(ebfh->unit);
	int32_t sdefsize = strlen(ebfh->sdef);
	int temp = 0;
	char pad = 60;
	int i,j;
	ebfh->ecode=0;
	if(ebfh->rank > EBF_RANK_MAXSIZE)
	{
		ebfh->ecode=311;
		ebfh->rank=EBF_RANK_MAXSIZE;
		printf("%s \n","Ebf Error: rank too large");
	}

	nwrit *= fwrite(&(ebfh->sig[0]), 8, 1, fp);
	nwrit *= fwrite(&(ebfh->version[0]), 4, 1, fp);

	unused4[0]=ebfh->endiantest;
	unused4[1]=ebfh->headersize;
	unused4[2]=namesize;
	unused4[3]=ebfh->datatype;
	unused4[4]=ebfh->datasize;
	unused4[5]=ebfh->rank;
	unused4[6]=unitsize;
	unused4[7]=sdefsize;
	if(ebfh->flag_swap ==1)
		ebfutils_SwapEndian(unused4,4,8);
	nwrit *= fwrite(unused4, 4*8, 1, fp);
	nwrit *= fwrite(&(ebfh->flags[0]), 4, 1, fp);


	i=0;
	unused8[i]=ebfh->capacity_;i++;
	for(j=0;j<ebfh->rank;++j)
	{
		unused8[i]=ebfh->dim[j];
		i++;
	}
	if(ebfh->flag_swap ==1)
		ebfutils_SwapEndian(unused8,8,i);
	nwrit *= fwrite(unused8, 8*i, 1, fp);


	if (namesize > 0)
		nwrit *= fwrite(ebfh->name, namesize, 1, fp);
	if (unitsize > 0)
	{
		nwrit *= fwrite(ebfh->unit, unitsize, 1, fp);
	}
	if (sdefsize > 0)
	{
		nwrit *= fwrite(ebfh->sdef, sdefsize, 1, fp);
	}

	temp = ebfh->headersize - (ftell(fp) - ebfh->headerpos);
	if (temp < 0)
	{
		ebfh->ecode=312;
		printf("%s \n","Ebf Error from EbfHeader::writeVer11()- exceeding headersize");
	}
	else if (temp > 0)
	{
		for (i = 0; i < temp; ++i)
			nwrit *= fwrite(&pad, 1, 1, fp);
	}

	if (nwrit != 1)
	{
		ebfh->ecode=306;
		printf("%s \n","Ebf Error from EbfHeader::writeVer11()- fwrite()");
	}

}


void EbfHeader_write(EbfHeader* ebfh, FILE* fp)
{
	ebfh->ecode=0;
	ebfh->headerpos = ftell(fp);
	if((ebfh->version[0]==1)&&(ebfh->version[1]==1))
		EbfHeader_writeVer11(ebfh, fp);
	else
	{
		printf("%s \n","Ebf Error from EbfHeader::write()- unrecognized header version");
		ebfh->ecode=313;
	}
	ebfh->datapos = ftell(fp);

}


static void EbfHeader_rename(EbfHeader* ebfh, FILE* fp,const char *name)
{
	int extrasize;
	ebfh->ecode=0;
	if((ebfh->version[0]==1)&&(ebfh->version[1]==1))
	{
		extrasize = ebfh->headersize
				- (56 + strlen(name) + strlen(ebfh->unit) + strlen(ebfh->sdef)
						+ 8 * ebfh->rank);

		if (extrasize >= 2)
		{
			if(ebfstrcpy(ebfh->name,EBF_KEYNAME_MAXSIZE, name)==0)
				EbfHeader_write(ebfh, fp);
			else
				ebfh->ecode = 307;
		}
		else
		{
			ebfh->ecode = 316;
		}
	}
	else
	{
		printf("%s %d%s%d \n","Ebf Error, for rename incompatible header version -> ",ebfh->version[0],".",ebfh->version[1]);
		ebfh->ecode = 317;
	}
}

static void EbfHeader_readVer10(EbfHeader* ebfh, FILE* fp)
{
	int i = 0;
	size_t j = 0;
	int64_t unused8[8];
	int fcheck = 0;
	char tempc1[101];
	memset(tempc1,0,sizeof(tempc1));
	ebfh->ecode=0;

	ebfh->flag_swap = 0;
	ebfh->headerpos = ftell(fp);
	ebfh->headersize = 256;
	fcheck = fread(&(ebfh->sig[0]), 6, 1, fp);
	fcheck = fread(tempc1, 100, 1, fp);

	ebfstrtrim(tempc1,sizeof(tempc1));
	if(ebfstrcpy(ebfh->name, EBF_KEYNAME_MAXSIZE,tempc1)!=0)
		ebfh->ecode=7;

	for(j=0;j<(strlen(ebfh->name));++j)
		ebfh->name[j]=tolower(ebfh->name[j]);


	fcheck = fread(&(ebfh->version[0]), 2, 1, fp);
	fcheck = fread(tempc1, 36, 1, fp);
	fcheck = fread(&(ebfh->endiantest), 4, 1, fp);
	fcheck = fread(&(ebfh->datatype), 4, 1, fp);
	fcheck = fread(&(ebfh->datasize), 4, 1, fp);
	fcheck = fread(&(ebfh->rank), 4, 1, fp);
	fcheck = fread(tempc1, 32, 1, fp);
	if (ebfh->endiantest != 256)
	{
		ebfh->flag_swap = 1;
		ebfutils_SwapEndian(&(ebfh->endiantest), 4, 1);
		ebfutils_SwapEndian(&(ebfh->datatype), 4, 1);
		ebfutils_SwapEndian(&(ebfh->datasize), 4, 1);
		ebfutils_SwapEndian(&(ebfh->rank), 4, 1);
	}
	if (ebfh->endiantest != 256)
	{
		printf("%s \n","Ebf Error from EbfHeader::ReadVer10()-cannot read header");
		ebfh->ecode=314;
	}
	else
	{
		if((ebfh->rank > EBF_RANK_MAXSIZE)||(ebfh->rank < 1))
		{
			ebfh->rank=EBF_RANK_MAXSIZE;
			ebfh->ecode=314;
		}

		if(fread(&unused8[0], 64, 1, fp) !=1)
			ebfh->ecode=306;
		for (i = 0; i < ebfh->rank; ++i)
			ebfh->dim[i] = unused8[i];

		if (ebfh->flag_swap == 1)
			ebfutils_SwapEndian(&(ebfh->dim[0]), 8, ebfh->rank);



		ebfh->datapos = ftell(fp);
		ebfh->unit[0] = 0;
		ebfh->sdef[0] = 0;
		ebfh->capacity_ = EbfHeader_elements(ebfh)*ebfh->datasize;

		ebfh->version[0] = 1;
		ebfh->version[1] = 0;
		ebfh->version[2] = 0;
		ebfh->version[3] = 0;
		ebfh->flags[0] = 0;
		ebfh->flags[1] = 0;
		ebfh->flags[2] = 0;
		ebfh->flags[3] = 0;
		ebfh->sig[0] = -118;
		ebfh->sig[1] = 69;
		ebfh->sig[2] = 66;
		ebfh->sig[3] = 70;
		ebfh->sig[4] = -82;
		ebfh->sig[5] = 43;
		ebfh->sig[6] = -81;
		ebfh->sig[7] = 10;
		ebfh->endiantest = 1684234849;


		if (fcheck != 1)
		{
			ebfh->ecode=306;
			printf("%s \n","Ebf Error from EbfHeader::ReadVer10()- fread()");
		}

		if ((ebfh->headerpos + ebfh->headersize) != ebfh->datapos)
		{
			ebfh->ecode=315;
			printf("%s \n","Ebf Error from EbfHeader::ReadVer10()-data position not correct in header");
		}
	}
}

static void EbfHeader_readVer11(EbfHeader* ebfh, FILE* fp)
{
	int32_t unused4[8];
	int namesize = 0;
	int unitsize = 0;
	int sdefsize = 0;
	int fcheck = 1;
	size_t i;
	ebfh->ecode=0;
	ebfh->flag_swap = 0;
	ebfh->headerpos = ftell(fp);
	fcheck *= fread(&(ebfh->sig[0]), 8, 1, fp);
	fcheck *= fread(&(ebfh->version[0]), 4, 1, fp);
	fcheck *= fread(&unused4[0], 4*8, 1, fp);
	fcheck *= fread(&(ebfh->flags[0]), 4, 1, fp);
	fcheck *= fread(&(ebfh->capacity_), 8, 1, fp);
	ebfh->endiantest=unused4[0];


	if (ebfh->endiantest != 1684234849)
	{
		ebfh->flag_swap = 1;
		ebfutils_SwapEndian(&(ebfh->endiantest), 4, 1);
	}
	if (ebfh->endiantest != 1684234849)
	{
		printf("%s \n","Ebf Error from EbfHeader::ReadVer11()-cannot read header");
		ebfh->ecode=314;
	}
	else
	{
		if (ebfh->flag_swap == 1)
		{
			ebfutils_SwapEndian(&unused4[0], 4, 8);
			ebfutils_SwapEndian(&(ebfh->capacity_), 8, 1);
		}
		ebfh->endiantest=unused4[0];
		ebfh->headersize=unused4[1];
		namesize=unused4[2];
		ebfh->datatype=unused4[3];
		ebfh->datasize=unused4[4];
		ebfh->rank=unused4[5];
		unitsize=unused4[6];
		sdefsize=unused4[7];

		if((ebfh->rank > EBF_RANK_MAXSIZE)||(ebfh->rank < 0))
		{
			ebfh->ecode=311;
		}
		if( ((int)(EBF_KEYNAME_MAXSIZE) <= (namesize)) || ((int)(EBF_KEYNAME_MAXSIZE) <= (unitsize)) || ((int)(EBF_SDEF_MAXSIZE) <= (sdefsize)) )
		{
				printf("%s \n","char size is too small");
				ebfh->ecode=311;
		}
		if(ebfh->ecode==0)
		{
			if(ebfh->rank>0)
			{
				fcheck *= fread(&(ebfh->dim[0]), 8 * ebfh->rank, 1, fp);
				if (ebfh->flag_swap == 1)
					ebfutils_SwapEndian(&(ebfh->dim[0]), 8, ebfh->rank);
			}
			else
			{
				ebfh->rank=1;
				ebfh->dim[0]=1;
			}

			if (namesize > 0)
				fcheck *= fread(ebfh->name, namesize, 1, fp);
			ebfh->name[namesize] = 0;
			for(i=0;i<(strlen(ebfh->name));++i)
				ebfh->name[i]=tolower(ebfh->name[i]);

			if (unitsize > 0)
			{
				fcheck *= fread(ebfh->unit, unitsize, 1, fp);
			}
			ebfh->unit[unitsize] = 0;

			if (sdefsize > 0)
			{
				fcheck *= fread(ebfh->sdef, sdefsize, 1, fp);
			}
			ebfh->sdef[sdefsize] = 0;
		}

		ebfh->extrasize = ebfh->headerpos + ebfh->headersize - ftell(fp);
		fseek(fp, (ebfh->extrasize), SEEK_CUR);

		ebfh->datapos = ftell(fp);
		if (fcheck != 1)
		{
			ebfh->ecode=306;
			printf("%s \n","Ebf Error from EbfHeader::ReadVer11()- fread()");
		}

	}

}

void EbfHeader_read(EbfHeader* ebfh, FILE* fp)
{
	char tempc1[12];
	char sig[8];
	ebfh->ecode=0;
	if (fread(tempc1, 12, 1, fp) == 1)
	{
		fseek(fp, -12, SEEK_CUR);
		/*
		//		if (strncmp("EBF>>>", tempc1, 6) == 0)
		//		{
		//			EbfHeader_readVer10(ebfh, fp);
		//		}
		//		else
		//		{
		//			if ((tempc1[8] == 1) && (tempc1[9] == 1))
		//				EbfHeader_readVer11(ebfh, fp);
		//			else
		//			{
		//				printf("%s \n","Ebf Error from EbfHeader::Read() unrecognized header version");
		//				ebfh->ecode=313;
		//			}
		//		}
		 */
		sig[0] = -118;
		sig[1] = 69;
		sig[2] = 66;
		sig[3] = 70;
		sig[4] = -82;
		sig[5] = 43;
		sig[6] = -81;
		sig[7] = 10;
		if (strncmp(sig, tempc1, 6) == 0)
		{
			if ((tempc1[8] == 1) && (tempc1[9] == 1))
				EbfHeader_readVer11(ebfh, fp);
			else
			{
				printf("%s \n","Ebf Error from EbfHeader::read() unrecognized header version");
				ebfh->ecode=313;
			}
		}
		else
		{
/*			if (strncmp("EBF>>>", tempc1, 6) == 0) */
			if (strncmp("EBF", tempc1, 3) == 0)
				EbfHeader_readVer10(ebfh, fp);
			else
			{
				printf("%s \n","Ebf Error from EbfHeader::read() unrecognized header version");
				ebfh->ecode=313;
			}

		}


	}
	else
	{
		ebfh->ecode=306;
		printf("%s \n","error reading ebf header");
	}

	if(ebfh->ecode != 0)
	{
		ebfh->rank = 1;
		ebfh->dim[0] = 0;
		ebfh->name[0] = 0;
		ebfh->unit[0] = 0;
		ebfh->sdef[0] = 0;
		ebfh->capacity_ = 0;
		ebfh->datatype = 0;
		ebfh->datasize = 0;
	}
	ebfh->datapos_end = ebfh->datapos + EbfHeader_elements(ebfh)*ebfh->datasize;
}

/*---------------------------------------------------------------*/

/**
 * Usagae internal: a structure to store data name keys and its locations
 */
typedef struct
{
	char key[EBF_KEYNAME_MAXSIZE];
	int64_t value;
} EbfKeyVal;


typedef struct
{
    char version[8];
    int64_t endiantest;
    int64_t headersize;
    int64_t typesize;
    int64_t hash_algo;
    int64_t count;
    int64_t htcapacity;
    int64_t nodesize;
    int64_t nodepos;
    int64_t nodecapacity;
    int64_t keypos;
    int64_t keycapacity;
    int64_t keyposcur;
} EbfFileHTHeader;

typedef struct
{
	int64_t keyloc;
    int64_t keysize;
    int64_t value;
    int64_t next;
    int64_t tnext;
} EbfFileHTNode;

/* void Ebf_FileHT_init_header(EbfFileHTHeader *header,int nodesize1,int capacity); */

typedef  struct
{
	int64_t hinfo[5];
	EbfFileHTHeader header;
	int flagswap,ecode,hecode;
	int64_t dpos_hinfo;
	EbfHeader ebfh_hinfo;
	FILE* fp;
} EbfFileHT;

/*----------------------------------------------------------------*/

/**
 * The djb2 hash function http://www.cse.yorku.ca/~oz/hash.html
 * This algorithm (k=33) was first reported by dan bernstein many years ago in comp.lang.c
 * @param mystr string to hash
 * @param hash1 previous hash if any to update or else uses value 5381
 * @return has value and 64 bit signed integer
 */
static int64_t EbfCkHash(const char *mystr, int64_t hash1)
{
	int64_t hash;
	size_t i, mylen;
	if (hash1 == 0)
	{
		hash = 5381;
	}
	else
	{
		hash = hash1;
	}
	mylen = strlen(mystr);
	for (i = 0; i < mylen; ++i)
		hash = hash * 33 + mystr[i];
	return hash;
}



static int64_t EbfLtHash(const char* mystr,int32_t capacity)
{
    uint64_t euhash=5381;
    size_t i,mylen;
	mylen=strlen(mystr);
	for(i=0;i<mylen;++i)
		euhash = euhash *33 + (uint64_t)(mystr[i]);
    return (euhash%capacity);

/*
    int64_t ehash=5381;
    size_t i,mylen;
	mylen=strlen(mystr);
	int64_t xmax=9223372036854775807;
	for(i=0;i<mylen;++i)
		ehash = ehash * 33 + (int64_t)(mystr[i]);
    ehash=ehash%capacity;
    if (ehash<0)
    {
    	ehash=(ehash+2*(xmax%capacity)+2)%capacity;
        if (ehash<0)
        {
        	ehash=ehash+capacity;
        }
    }
    return ehash;
 */
}

static void EbfFileHTHeader_Init(EbfFileHTHeader *header,int nodesize1,int capacity)
{
	int i;
    for(i=0;i<8;++i)
    	header->version[i]=0;
    header->version[0]=1;
    header->endiantest=1684234849;
    header->headersize=sizeof(header[0]);
    header->typesize=sizeof(header->count);
    header->hash_algo=1;
    header->count=0;
    header->htcapacity=capacity+capacity/2;
    header->nodesize=nodesize1;
    header->nodepos=header->headersize+header->typesize*header->htcapacity;
    header->nodecapacity=capacity;
    header->keypos=header->nodepos+header->nodesize*header->nodecapacity;
    header->keycapacity=20*header->nodecapacity;
    header->keyposcur=0;
}

/*
static void EbfTable_write_key(EbfFileHT *fileht, int64_t loc,const char*  key)
{
	if(fileht->ecode ==0)
	{
		if(strlen(key)>=EBF_KEYNAME_MAXSIZE)
		{
			printf("EBF Error: write_key() keysize too large ");
			fileht->ecode=207;
		}
		else
		{
			fseek(fileht->fp,fileht->hinfo[2]+fileht->header.keypos+loc,SEEK_SET);
			if(fwrite(key,strlen(key), 1,fileht->fp) != 1)
				fileht->ecode=206;
		}
	}
}
*/

static void EbfFileHT_read_key(EbfFileHT *fileht, int64_t loc, int64_t keysize,char* key)
{
	size_t i;
	if(fileht->ecode ==0)
	{
		char temp[EBF_KEYNAME_MAXSIZE];
		if((keysize+1)>(int)(EBF_KEYNAME_MAXSIZE))
		{
			printf("%s %d \n","EBF Error: FileHT_read_key() keysize too large ",(int)keysize);
			fileht->ecode=207;
		}
		else
		{
			fseek(fileht->fp,fileht->hinfo[2]+fileht->header.keypos+loc,SEEK_SET);
			if(fread(temp,keysize, 1,fileht->fp) == 1)
			{
				strncpy(key,temp,keysize);
				key[keysize]=0;
				for(i=0;i<strlen(key);++i)
					key[i]=tolower(key[i]);
			}
			else
			{
				fileht->ecode=206;
			}
		}
	}
}

static void EbfFileHT_write_header(EbfFileHT *fileht)
{
	if(fileht->ecode ==0)
	{
		fseek(fileht->fp,fileht->hinfo[2],SEEK_SET);
		if(fwrite(&(fileht->header.version), 8, 1, fileht->fp) != 1)
			fileht->ecode=206;
		if (fileht->flagswap==1)
		{
			ebfutils_SwapEndian(&(fileht->header.endiantest), sizeof(fileht->header.endiantest),12);
			if(fwrite(&(fileht->header.endiantest), sizeof(fileht->header.endiantest)*12, 1, fileht->fp) !=1)
				fileht->ecode=206;
			ebfutils_SwapEndian(&(fileht->header.endiantest), sizeof(fileht->header.endiantest),12);
		}
		else
		{
			if(fwrite(&(fileht->header.endiantest), sizeof(fileht->header.endiantest)*12, 1, fileht->fp) != 1)
				fileht->ecode=206;
		}
	}
}


static void EbfFileHT_write_header_gen(EbfFileHTHeader *header,FILE* fp,int flagswap,int *ecode)
{
		if(fwrite(&(header->version), 8, 1, fp) != 1)
			(*ecode)=206;
		if (flagswap==1)
		{
			ebfutils_SwapEndian(&(header->endiantest), sizeof(header->endiantest),12);
			if(fwrite(&(header->endiantest), sizeof(header->endiantest)*12, 1, fp) !=1)
				(*ecode)=206;
			ebfutils_SwapEndian(&(header->endiantest), sizeof(header->endiantest),12);
		}
		else
		{
			if(fwrite(&(header->endiantest), sizeof(header->endiantest)*12, 1, fp) != 1)
				(*ecode)=206;
		}
}


static void EbfFileHT_read_header(EbfFileHT *fileht)
{
	if(fileht->ecode ==0)
	{
		fileht->header.version[0]=0;
		fileht->header.endiantest=0;
		fileht->flagswap=0;
		fseek(fileht->fp,fileht->hinfo[2],SEEK_SET);
		if(fread(&(fileht->header.version), 8, 1, fileht->fp) != 1)
			fileht->ecode=206;
		if(fread(&(fileht->header.endiantest), sizeof(fileht->header.endiantest)*12, 1, fileht->fp) != 1)
			fileht->ecode=206;
		if(fileht->header.version[0] !=1)
			fileht->ecode=206;
		if(fileht->header.endiantest!=1684234849)
		{
			fileht->flagswap=1;
			ebfutils_SwapEndian(&(fileht->header.endiantest), sizeof(fileht->header.endiantest),12);
		}
/* retouch */
		if(fileht->header.endiantest!=1684234849)
			fileht->ecode=206;

	}
}
static void EbfFileHT_read_hvalue(EbfFileHT *fileht, int64_t loc,int64_t *hvalue)
{
	if(fileht->ecode ==0)
	{
		fseek(fileht->fp,fileht->hinfo[2]+fileht->header.headersize+fileht->header.typesize*loc,SEEK_SET);
		if(fread(hvalue,fileht->header.typesize,1,fileht->fp) != 1)
			fileht->ecode=206;
		if(fileht->flagswap==1)
			ebfutils_SwapEndian(hvalue, fileht->header.typesize,1);
	}
}
static void EbfFileHT_write_hvalue(EbfFileHT *fileht, int64_t loc,int64_t *hvalue)
{
	if(fileht->ecode ==0)
	{
		fseek(fileht->fp,fileht->hinfo[2]+fileht->header.headersize+fileht->header.typesize*loc,SEEK_SET);
		if(fileht->flagswap==1)
		{
			ebfutils_SwapEndian(hvalue, fileht->header.typesize,1);
			if(fwrite(hvalue,fileht->header.typesize,1,fileht->fp) != 1)
				fileht->ecode=206;
			ebfutils_SwapEndian(hvalue, fileht->header.typesize,1);
		}
		else
		{
			if( fwrite(hvalue,fileht->header.typesize,1,fileht->fp) != 1)
				fileht->ecode=206;
		}
	}
}

static void EbfFileHT_read_node(EbfFileHT *fileht, int64_t loc,EbfFileHTNode *nodep)
{
	if(fileht->ecode ==0)
	{
		fseek(fileht->fp,fileht->hinfo[2]+fileht->header.nodepos+fileht->header.nodesize*loc,SEEK_SET);
		if(fread(nodep,fileht->header.nodesize,1,fileht->fp) != 1)
			fileht->ecode=206;
		if(fileht->flagswap==1)
			ebfutils_SwapEndian(nodep, fileht->header.typesize,fileht->header.nodesize/fileht->header.typesize);
	}
}
static void EbfFileHT_write_node(EbfFileHT *fileht,int64_t loc,EbfFileHTNode *nodep)
{
	if(fileht->ecode ==0)
	{
		fseek(fileht->fp,fileht->hinfo[2]+fileht->header.nodepos+fileht->header.nodesize*loc,SEEK_SET);
		if(fileht->flagswap==1)
		{
			ebfutils_SwapEndian(nodep, fileht->header.typesize,fileht->header.nodesize/fileht->header.typesize);
			if(fwrite(nodep,fileht->header.nodesize,1,fileht->fp) != 1)
				fileht->ecode=206;
			ebfutils_SwapEndian(nodep, fileht->header.typesize,fileht->header.nodesize/fileht->header.typesize);
		}
		else
		{
			if(fwrite(nodep,fileht->header.nodesize,1,fileht->fp) != 1)
				fileht->ecode=206;
		}
	}
}

static void EbfFileHT_write_hinfo(EbfFileHT *fileht)
{
/*	if(fileht->ecode ==0)
	{*/
	fseek(fileht->fp,fileht->dpos_hinfo,SEEK_SET);
	if (fileht->ebfh_hinfo.flag_swap == 1)
	{
		ebfutils_SwapEndian(&(fileht->hinfo[0]),sizeof(fileht->hinfo[0]),EbfHeader_elements(&(fileht->ebfh_hinfo)));
		if(fwrite(&(fileht->hinfo[0]),sizeof(fileht->hinfo[0])*EbfHeader_elements(&(fileht->ebfh_hinfo)),1,fileht->fp) != 1)
			fileht->ecode=206;
		ebfutils_SwapEndian(&(fileht->hinfo[0]),sizeof(fileht->hinfo[0]),EbfHeader_elements(&(fileht->ebfh_hinfo)));
	}
	else
	{
		if(fwrite(&(fileht->hinfo[0]),sizeof(fileht->hinfo[0])*EbfHeader_elements(&(fileht->ebfh_hinfo)),1,fileht->fp) != 1)
			fileht->ecode=206;
	}
/*	} */
}


static void EbfFileHT_read_hinfo(EbfFileHT *fileht)
{
	if(fileht->ecode ==0)
	{
		fseek(fileht->fp,fileht->dpos_hinfo,SEEK_SET);
		if(fread(&(fileht->hinfo[0]),sizeof(fileht->hinfo[0])*EbfHeader_elements(&(fileht->ebfh_hinfo)),1,fileht->fp) != 1)
			fileht->ecode=206;
		if (fileht->ebfh_hinfo.flag_swap == 1)
			ebfutils_SwapEndian(&(fileht->hinfo[0]),sizeof(fileht->hinfo[0]),EbfHeader_elements(&(fileht->ebfh_hinfo)));
	}
}


static void EbfFileHT_close(EbfFileHT* fileht)
{
	if(fileht->fp != NULL)
	{
		fclose(fileht->fp);
		fileht->fp=NULL;
	}
}


/*static void EbfFileHT_create_htable(EbfFileHT *fileht, int64_t capacity) */
static void EbfFileHT_create_htable(const char* filename, int64_t capacity)
{
	int i,nsize,nwrite;
	int64_t dims[8];
	EbfFileHTNode node;
	EbfFileHTNode* nodearr=NULL;
	char* keyarr=NULL;
	int64_t* htarr=NULL;
	EbfHeader ebfh,ebfh_hinfo;
	EbfFileHTHeader header;
	int64_t hinfo[10];
	int64_t dpos_hinfo;
	int ecode=0;
	FILE* fp=NULL;

	if (capacity>1000000)
	{
		printf("%s \n","capacity>10^6");
	}

	if(ecode ==0)
	{
		fp=fopen(filename,"rb+");
		fseek(fp,0,SEEK_SET);
		EbfHeader_read(&(ebfh_hinfo), fp);
		dpos_hinfo = ftell(fp);
		ecode = ebfh_hinfo.ecode;
		if ((strcmp(ebfh_hinfo.name, "/.ebf/info") != 0)||(ebfh_hinfo.datatype != 3)||(EbfHeader_elements(&(ebfh_hinfo)) > 10))
		{
				ecode = 216;
				printf("%s \n","EBF Error from EbfHT: unrecognized data type for hinfo or size");
		}

		if (0==ecode)
		{
			if(fread(&(hinfo[0]),sizeof(hinfo[0])*EbfHeader_elements(&(ebfh_hinfo)),1,fp) != 1)
				ecode=206;
			if (ebfh_hinfo.flag_swap == 1)
				ebfutils_SwapEndian(&(hinfo[0]),sizeof(hinfo[0]),EbfHeader_elements(&(ebfh_hinfo)));
		}
		if (0==ecode)
		{
			fseek(fp,0,SEEK_END);
			hinfo[1]=ftell(fp);
			EbfFileHTHeader_Init(&(header),sizeof(node),capacity);
			dims[0]=(sizeof(header)+header.typesize*header.htcapacity+header.nodesize*header.nodecapacity+header.keycapacity);
			EbfHeader_set(&ebfh,"/.ebf/htable",1,1,"","",1,dims);

			/* added*/
			ebfh.flag_swap=ebfh_hinfo.flag_swap;

			EbfHeader_write(&ebfh,fp);
			hinfo[2]=ftell(fp);
			hinfo[3]=1;

			EbfFileHT_write_header_gen(&header,fp,ebfh.flag_swap,&ecode);

			htarr=(int64_t* )malloc(sizeof(int64_t)*header.htcapacity);
			if((htarr != NULL)&&(header.typesize==sizeof(int64_t)))
			{
				for(i=0;i<header.htcapacity;i++)
					htarr[i]=0;

				/* added*/
				if(ebfh.flag_swap==1)
					ebfutils_SwapEndian(htarr, header.typesize,header.htcapacity);

				if(fwrite(htarr,header.typesize*header.htcapacity,1,fp)!=1)
					ecode=206;
				free(htarr);
				htarr=NULL;
			}
			else
				ecode=204;

			node.keyloc=0;
			node.keysize=0;
			node.value=-1;
			node.next=-1;
			node.tnext=-1;
			nodearr=(EbfFileHTNode* )malloc(sizeof(EbfFileHTNode)*header.nodecapacity);
			if((nodearr != NULL)&&(sizeof(EbfFileHTNode)==header.nodesize))
			{
				for(i=0;i<header.nodecapacity;i++)
					nodearr[i]=node;
				/* added*/
				if(ebfh.flag_swap==1)
					ebfutils_SwapEndian(nodearr, header.typesize,header.nodecapacity*header.nodesize/header.typesize);

				if(fwrite(nodearr,header.nodesize*header.nodecapacity,1,fp)!= 1)
					ecode=206;
				free(nodearr);
				nodearr=NULL;
			}
			else
				ecode=204;

			keyarr=(char* )malloc(header.keycapacity);
			if(keyarr != NULL)
			{
				for(i=0;i<header.keycapacity;i++)
					keyarr[i]=0;
				if(fwrite(keyarr,header.keycapacity,1,fp) != 1)
					ecode=206;
				free(keyarr);
				keyarr=NULL;
			}
			else
				ecode=204;

			nsize=dims[0]-(ftell(fp)-hinfo[2]);
			if(nsize != 0)
			{
				nwrite=fwrite(&nwrite,nsize,1,fp);
				ecode=211;
			}
			if(0==ecode)
			{
				fseek(fp,dpos_hinfo,SEEK_SET);
				if (ebfh_hinfo.flag_swap == 1)
				{
					ebfutils_SwapEndian(&(hinfo[0]),sizeof(hinfo[0]),EbfHeader_elements(&(ebfh_hinfo)));
					if(fwrite(&(hinfo[0]),sizeof(hinfo[0])*EbfHeader_elements(&(ebfh_hinfo)),1,fp) != 1)
						ecode=206;
					ebfutils_SwapEndian(&(hinfo[0]),sizeof(hinfo[0]),EbfHeader_elements(&(ebfh_hinfo)));
				}
				else
				{
					if(fwrite(&(hinfo[0]),sizeof(hinfo[0])*EbfHeader_elements(&(ebfh_hinfo)),1,fp) != 1)
						ecode=206;
				}
			}

		}

		if(fp!=NULL)
		{
			fclose(fp);
			fp=NULL;
		}

	}

}



static void EbfFileHT_addKey(EbfFileHT *fileht,const char* key, int64_t value)
{
	EbfFileHTNode node;
	int64_t hindex,hvalue,loc;
	if(fileht->ecode == 0)
	{
		fileht->header.count++;
		if((fileht->header.count+1) > fileht->header.nodecapacity)
		{
			printf("%s \n","Ebf error: not enough space for more keys");
			fileht->ecode=211;
		}
		if((fileht->header.keyposcur+(int)strlen(key)+1) > fileht->header.keycapacity)
		{
			printf("%s \n","Ebf error: not enough space for more keys");
			fileht->ecode=211;
		}
		if(fileht->ecode == 0)
		{

		node.keyloc=fileht->header.keyposcur;
		node.keysize=strlen(key);
		node.value=value;
		node.next=-1;
		node.tnext=-1;



		fseek(fileht->fp,fileht->hinfo[2]+fileht->header.keypos+node.keyloc,SEEK_SET);
		if(fwrite(key,strlen(key), 1,fileht->fp) !=1 )
			fileht->ecode=206;
		EbfFileHT_write_node(fileht, fileht->header.count,&node);

		fileht->header.keyposcur+=node.keysize;
		EbfFileHT_write_header(fileht);

		hindex=EbfLtHash(key,fileht->header.htcapacity);
		EbfFileHT_read_hvalue(fileht, hindex,&hvalue);
		if(fileht->ecode == 0)
		{
			if (hvalue != 0)
			{
				loc=hvalue;
				EbfFileHT_read_node(fileht, loc,&node);
				while(node.next != -1)
				{
					loc=node.next;
					EbfFileHT_read_node(fileht, loc,&node);
					if(fileht->ecode != 0)
						break;
				}
				node.next=fileht->header.count;
				EbfFileHT_write_node(fileht, loc,&node);
			}
			else
			{
				EbfFileHT_write_hvalue(fileht, hindex,&(fileht->header.count));
			}
		}
		}
	}
}

/* only routine that takes fileht and initializes fileht->ecode */
/*  others only work if fileht->ecode is 0*/
static void EbfFileHT_setup(EbfFileHT *fileht, const char* filename, const char* mode)
{
	int i;
	fileht->fp = NULL;
	fileht->flagswap = 0;
	for (i = 0; i < 10; ++i)
		fileht->hinfo[i] = 0;
	fileht->fp = fopen(filename, mode);
	fileht->hecode = 0;
	fileht->ecode = 0;

	if (fileht->fp != NULL)
	{

		EbfHeader_read(&(fileht->ebfh_hinfo), fileht->fp);
		fileht->dpos_hinfo = ftell(fileht->fp);
		fileht->ecode = fileht->ebfh_hinfo.ecode;
/*		fileht->hecode = fileht->ebfh_hinfo.ecode; */

		if ((strcmp(fileht->ebfh_hinfo.name, "/.ebf/info") != 0) ||
			(fileht->ebfh_hinfo.datatype != 3) ||
			(EbfHeader_elements(&(fileht->ebfh_hinfo)) > 10)
			|| (EbfHeader_elements(&(fileht->ebfh_hinfo)) < 4))
		{
			fileht->ecode = 216;
			fileht->hecode = 216;
		}


		if (fileht->ecode == 0)
		{
			if (fread(
					&(fileht->hinfo[0]),
					sizeof(fileht->hinfo[0])
							* EbfHeader_elements(&(fileht->ebfh_hinfo)), 1,
					fileht->fp) != 1)
			{
				printf("%s \n", "EBF error in fread");
				fileht->ecode = 206;
				fileht->hecode = 206;
			}
			else
			{
				/* this is where hecode differs from ecode*/
				/* .e.g %/.ebf/htable might not exist but /.ebf/info exists */
				if (fileht->ebfh_hinfo.flag_swap)
					ebfutils_SwapEndian(&(fileht->hinfo[0]),
							sizeof(fileht->hinfo[0]),
							EbfHeader_elements(&(fileht->ebfh_hinfo)));
/*					EbfHeader_read(&(fileht->ebfh_htable), fileht->fp); */
				if((fileht->hinfo[1]>0)&&(fileht->hinfo[2]>0)&&(fileht->hinfo[3]==1))
				{
					EbfFileHT_read_header(fileht);
				}
				else
				{
					fileht->hecode = 211;
				}
			}
		}
	}
	else
	{
		fileht->ecode = 205;
		fileht->hecode = 205;
	}
}


static void EbfTable_getKeyValHT(const char* filename, EbfKeyVal **keyvalp,int *kv_size,int* kv_capacity, int* ecode)
{
	int i,nread;
	int64_t hvalue,loc;
	int64_t* htable=NULL;
	char* keyarr=NULL;
	EbfFileHTNode* nodearr=NULL;
	int ncapacity,nsize;
	EbfFileHTNode node;
	EbfKeyVal* temp_p;
	char key1[EBF_KEYNAME_MAXSIZE];
	EbfFileHT fileht1;
	EbfFileHT* fileht;
	fileht=&fileht1;
	fileht->ecode=0;
	EbfFileHT_setup(fileht,filename,"rb");
	ncapacity=0;
	nsize=0;
	nread=1;

	if (((*keyvalp) != NULL)||(fileht->header.typesize != 8)||(fileht->ecode!=0))
	{
		printf("%s \n", "EBF error input arrays should be NULL or");
		printf("%s \n","typesize problem in EbfHT or setup problem");
		fileht->ecode=218;
	}
	else
	{
		htable=NULL;
		fseek(fileht->fp,fileht->hinfo[2]+fileht->header.headersize,SEEK_SET);
		htable=(int64_t *)malloc(sizeof(int64_t)*fileht->header.htcapacity);
		if((htable != NULL)&&(fileht->header.typesize==sizeof(int64_t)))
		{
			nread*=fread(htable,fileht->header.typesize*fileht->header.htcapacity,1,fileht->fp);
			if(fileht->flagswap==1)
				ebfutils_SwapEndian(htable, fileht->header.typesize*fileht->header.htcapacity,1);
		}
		else
			nread=0;
		keyarr=NULL;
		fseek(fileht->fp,fileht->hinfo[2]+fileht->header.keypos,SEEK_SET);
		keyarr=(char *)malloc(fileht->header.keycapacity);
		if(keyarr != NULL)
		{
			nread*=fread(keyarr,fileht->header.keycapacity,1,fileht->fp);
		}
		else
			nread=0;

		nodearr=NULL;
		fseek(fileht->fp,fileht->hinfo[2]+fileht->header.nodepos,SEEK_SET);
		nodearr=(EbfFileHTNode *)malloc(sizeof(EbfFileHTNode)*fileht->header.nodecapacity);
		if((nodearr != NULL)&&(sizeof(EbfFileHTNode)==fileht->header.nodesize))
		{
			nread*=fread(nodearr,fileht->header.nodesize*fileht->header.nodecapacity,1,fileht->fp);
			if((fileht->flagswap==1)&&(nread==1))
				ebfutils_SwapEndian(nodearr, fileht->header.typesize,fileht->header.nodesize*fileht->header.nodecapacity/fileht->header.typesize);
		}
		else
			nread=0;

		if(nread != 1)
			fileht->ecode=206;

		if(fileht->ecode == 0)
		{

			for(i=0;i<fileht->header.htcapacity;i++)
			{
				hvalue=htable[i];
				if (hvalue != 0)
				{
					loc=hvalue;
					do
					{
						node=nodearr[loc];
						if(node.keysize >= (int)(EBF_KEYNAME_MAXSIZE))
						{
							printf("%s \n","node.keysize problem in EbfHT");
							fileht->ecode=207;
							key1[0]=0;
						}
						else
						{
							strncpy(key1,&keyarr[node.keyloc],node.keysize);
							key1[node.keysize]=0;
						}
						if (nsize==ncapacity)
						{
							ncapacity*=2;
							if(ncapacity==0)
								ncapacity=16;
							temp_p=(EbfKeyVal* )realloc(*keyvalp, ncapacity*sizeof(EbfKeyVal));
							if(temp_p != NULL)
								*keyvalp=temp_p;
							else
								fileht->ecode=204;

						}
						strcpy((*keyvalp)[nsize].key,key1);

						(*keyvalp)[nsize].value=node.value;
						nsize++;
						loc=node.next;
						if(nsize > fileht->header.nodecapacity)
						{
							fileht->ecode=215;
						}
						if(fileht->ecode != 0)
							break;

					}
					while(node.next != -1);
				}
			}
			if(htable!=NULL)
			{
				free(htable);
				htable=NULL;
			}
			if(keyarr!=NULL)
			{
				free(keyarr);
				keyarr=NULL;
			}
			if(nodearr!=NULL)
			{
				free(nodearr);
				nodearr=NULL;
			}
		}
	}

	EbfFileHT_close(fileht);
	if(fileht->ecode != 0)
		nsize=0;
	(*kv_size)=nsize;
	(*kv_capacity)=ncapacity;
}


/* if fileht->ecode is zero means hashtable exists Read from there
 * else try iterative.
*/
static int64_t EbfFileHT_getfp(EbfFileHT* fileht, const char* key)
{
	int64_t loc,hindex,iter,temploc;
	EbfFileHTNode node;
	EbfHeader ebfh;
	char key1[EBF_KEYNAME_MAXSIZE];
	int64_t value=-1;

	if(fileht->ecode==0)
	{
		hindex=EbfLtHash(key,fileht->header.htcapacity);
		EbfFileHT_read_hvalue(fileht, hindex,&loc);
		if (loc != 0)
		{
			EbfFileHT_read_node(fileht, loc,&node);
			if(fileht->ecode!=0)
			{
				printf("%s \n","Ebf error,in EbfFileHT_getfp() -> read_node()");
/*				printf("%s %"PRId64" %"PRId64" %"PRId64" %"PRId64" \n","filename Error2: ",hindex,loc,fileht->header.nodecapacity,fileht->header.htcapacity); */
			}
			EbfFileHT_read_key(fileht, node.keyloc,node.keysize,key1);
			iter=0;
			while((node.next != -1)&&(strcmp(key1,key)!=0))
			{
				iter++;
				EbfFileHT_read_node(fileht, node.next,&node);
				EbfFileHT_read_key(fileht, node.keyloc,node.keysize,key1);
                if (iter > fileht->header.nodecapacity)
                {
                    printf("%s \n","Ebf erro: ebfht_getfp too many iterations ");
                    fileht->ecode=209 ;
                }
				if(fileht->ecode != 0)
					break;
			}
			if((strcmp(key1,key)==0)&&(fileht->ecode==0))
				value=node.value;
		}
	}
	else
	{
		fseek(fileht->fp,0,SEEK_END);
		loc=ftell(fileht->fp);
		fseek(fileht->fp,0,SEEK_SET);
		while(ftell(fileht->fp)<loc)
		{
			temploc=ftell(fileht->fp);
			EbfHeader_read(&ebfh,fileht->fp);
			if(strcmp(ebfh.name,key)==0)
			{
				value=temploc;
				break;
			}
			if(ebfh.ecode!=0)
			{
				fileht->ecode=ebfh.ecode;
				value=-2;
				break;
			}
			fseek(fileht->fp,ebfh.capacity_,SEEK_CUR);
		}
	}

	return value;
}


int64_t EbfTable_Get(const char* filename, const char* keyin,int* ecode)
{
	char key[EBF_KEYNAME_MAXSIZE];
	size_t i;
	EbfFileHT fileht1;
	EbfFileHT* fileht;
	int64_t value=-1;
	fileht=&fileht1;
	fileht->ecode=0;
	*ecode=0;
	EbfFileHT_setup(fileht, filename,"rb");
	if(fileht->fp != NULL)
	{	  
		if(ebfstrcpy(key,EBF_KEYNAME_MAXSIZE,keyin)==0)
		{
			for(i=0;i<strlen(key);++i)
				key[i]=tolower(key[i]);
			value=EbfFileHT_getfp(fileht,key);
			if(value==-2)
				*ecode=211;
		}
		else
		{
			printf("%s \n","Ebf error, key size too large ");
			*ecode=207;
		}		
	}
	else
	{
		printf("%s %s \n","Error opening file: ",filename);
		*ecode=205;
	}

	EbfFileHT_close(fileht);
	return value;
}

static int64_t EbfTable_Get_WithHeader(const char* filename, const char* keyin,EbfHeader *ebfh)
{
	char key[EBF_KEYNAME_MAXSIZE];
	size_t i;
	EbfFileHT fileht1;
	EbfFileHT* fileht;
	int64_t value=-1;
	fileht=&fileht1;
	fileht->ecode=0;
	ebfh->ecode=0;
	EbfFileHT_setup(fileht, filename,"rb");
	if(fileht->fp != NULL)
	{
		if(ebfstrcpy(key,EBF_KEYNAME_MAXSIZE,keyin)==0)
		{
			for(i=0;i<strlen(key);++i)
				key[i]=tolower(key[i]);
			value=EbfFileHT_getfp(fileht,key);
			if(value==-2)
				ebfh->ecode=211;
			if(value>=0)
			{
				fseek(fileht->fp,value,SEEK_SET);
				EbfHeader_read(ebfh,fileht->fp);
				if(ebfh->ecode==0)
				{
					if(strcmp(ebfh->name,key)!=0)
					{
						value=-1;
					}
				}
				else
					value=-1;

			}
		}
		else
		{
			printf("%s \n","Ebf error, key size too large ");
			ebfh->ecode=207;
		}
	}
	else
	{
		printf("%s %s \n","Error opening file: ",filename);
		ebfh->ecode=205;
	}

	EbfFileHT_close(fileht);
	return value;
}


/**
 * Will expand the existing table. Copies current elements in table, then creates a new
 * table and puts the copied values in.
 * @param filename
 * @param capacity
 * @param ecode
 */

static void EbfTable_expand(EbfFileHT *fileht,const char* filename,const char* keyin,int* ecode)
{
	EbfKeyVal* keyval=NULL;
	EbfHeader ebfh;
	int nsize=0;
	int ncapacity=0;
	int nexti=-1;
	int i;
	int capacity,temp,factor;
	char mykey[EBF_KEYNAME_MAXSIZE];
	int64_t myvalue;
	/*
	EbfFileHT fileht1;
	EbfFileHT *fileht;
	fileht=&fileht1;
	fileht->ecode=0;
	 */
	capacity=fileht->header.nodecapacity;
	factor=(fileht->header.keycapacity/fileht->header.nodecapacity);
	EbfFileHT_close(fileht);

	do
	{
		nexti++;
		sprintf(mykey,"%s%d","/.tr/.ebf/htable.t",nexti);
		if(strlen(mykey)>=EBF_KEYNAME_MAXSIZE)
		{
			fileht->ecode=207;
			printf("%s \n","Ebf error strlen(s)>100");
			break;
		}
	}
	while(EbfTable_Get(filename,mykey,&(fileht->ecode))!=-1);

	if(fileht->ecode == 0)
	{
		EbfTable_getKeyValHT(filename,&keyval,&nsize,&ncapacity,&(fileht->ecode));
	}
	if(fileht->ecode == 0)
	{
		EbfFileHT_setup(fileht, filename,"rb+");
		if(fileht->ecode == 0)
		{

			fseek(fileht->fp,fileht->hinfo[1],SEEK_SET);
			EbfHeader_read(&ebfh,fileht->fp);
			fseek(fileht->fp,fileht->hinfo[1],SEEK_SET);
			EbfHeader_rename(&ebfh,fileht->fp,mykey);
			fileht->ecode=ebfh.ecode;


			if(fileht->ecode==0)
			{
				myvalue=fileht->hinfo[1];
				temp=strlen(keyin)+strlen(mykey);
				for(i=0;i<nsize;++i)
				{
					temp+=strlen(keyval[i].key);
				}
				capacity=capacity*2;
				if (capacity<=0)
					capacity=16;
				while ( (nsize+3) > capacity)
				{
					capacity=capacity*2;
				}
				while ((temp+1) > capacity*factor)
				{
					capacity=capacity*2;
				}

				EbfFileHT_close(fileht);
				EbfFileHT_create_htable(filename, capacity);
				EbfFileHT_setup(fileht, filename,"rb+");

				for(i=0;i<nsize;++i)
				{
					if (strcmp(keyval[i].key,"/.ebf/htable")==0)
					{
						keyval[i].value=fileht->hinfo[1];
						break;
					}
				}

				for(i=0;i<nsize;++i)
				{
					EbfFileHT_addKey(fileht, keyval[i].key,keyval[i].value);
				}
				EbfFileHT_addKey(fileht, mykey,myvalue);
			}
		}
		EbfFileHT_close(fileht);
	}

	if(keyval!=NULL)
	{
		free(keyval);
		keyval=NULL;
	}
	EbfFileHT_setup(fileht, filename,"rb+");
	*ecode=fileht->ecode;
}



int EbfTable_Put(const char* filename,const char* keyin, int64_t value)
{
	int64_t capacity;
	size_t i;
	int ecode;
	char mystr[EBF_KEYNAME_MAXSIZE + 100];
	char key[EBF_KEYNAME_MAXSIZE];
	EbfFileHT fileht1;
	EbfFileHT* fileht;
	fileht=&fileht1;
	ecode=0;
	if(ebfstrcpy(key,EBF_KEYNAME_MAXSIZE,keyin)!=0)
	{
		ecode=207;
	}
	else
	{
		for (i = 0; i < strlen(key); ++i)
			key[i] = tolower(key[i]);

		EbfFileHT_setup(fileht, filename,"rb+");
		if(fileht->fp != NULL)
		{
			if(0==fileht->ecode)
			{
				if(EbfFileHT_getfp(fileht,key) < 0)
				{
					capacity=fileht->header.nodecapacity;
					if (capacity<=0)
						capacity=16;
					while ( (fileht->header.count+2) > capacity)
						capacity=capacity*2;
					while (((int)(strlen(key))+fileht->header.keyposcur+1) > capacity*(fileht->header.keycapacity/fileht->header.nodecapacity))
						capacity=capacity*2;
					if (capacity!=fileht->header.nodecapacity)
					{
						EbfTable_expand(fileht,filename,key,&(fileht->ecode));
					}

					if(0==fileht->ecode)
					{
						EbfFileHT_addKey(fileht, key,value);
					}
				}
				else
				{
					fileht->ecode=215;
				}
			}
			if(0==fileht->hecode)
			{
#ifdef EBFLANGOPTCPP
				std::stringstream mystr1;
				mystr1<<"("<<key<<", "<<value<<")";
				if (mystr1.str().size() >= (EBF_KEYNAME_MAXSIZE + 100))
				{
					fileht->ecode=207;
					printf("%s \n","Ebf Error from Efile.close()- mystr max size is too small");
				}
				else
				{
					fileht->hinfo[0] = EbfCkHash(mystr1.str().c_str(), fileht->hinfo[0]);
					EbfFileHT_write_hinfo(fileht);
					mystr[0]=0;
				}
#else

#ifdef EBFLANGOPTIDL
				/* On systems where PRId64 is not available */
				ebf_prid64_idl(value,value_str,EBF_KEYNAME_MAXSIZE);
				sprintf(mystr, "%s%s%s%s%s", "(", key, ", ", value_str, ")");
#else
				sprintf(mystr, "%s%s%s%"PRId64"%s", "(", key, ", ", value, ")");
#endif
	/*			sprintf(mystr, "%s%s%s %lld%s", "(", key, ",", value, ") "); */
				if (strlen(mystr) >= (EBF_KEYNAME_MAXSIZE + 100))
				{
					fileht->ecode=207;
					printf("%s \n","Ebf Error from Efile.close()- mystr max size is too small");
				}
				else
				{
					fileht->hinfo[0] = EbfCkHash(mystr, fileht->hinfo[0]);
					EbfFileHT_write_hinfo(fileht);
				}
#endif

			}
			ecode=fileht->ecode;
		}
		else
		{
			printf("%s %s \n","Error opening file: ",filename);
			ecode=205;
		}
		EbfFileHT_close(fileht);
	}

return ecode;
}



int EbfTable_Init(const char* filename)
{
	int i;
	int64_t dims[8];
	EbfHeader ebfh;
	EbfFileHT fileht1;
	EbfFileHT* fileht;
	fileht=&fileht1;
	fileht->ecode=0;
	fileht->fp = fopen(filename, "wb");
	if (fileht->fp != NULL)
	{
		dims[0] = 5;
		EbfHeader_set(&ebfh, "/.ebf/info", 3, 8, "", "", 1, dims);
		for (i = 0; i < 5; ++i)
			fileht->hinfo[i] = 0;
		/* just making sure capacity_ is same as datasize, if in future _set does not do it*/
		ebfh.capacity_=EbfHeader_elements(&ebfh)*ebfh.datasize;
		EbfHeader_write(&ebfh, fileht->fp);
		fileht->ecode=ebfh.ecode;
		if(fwrite(fileht->hinfo,ebfh.capacity_, 1, fileht->fp) != 1)
			fileht->ecode=206;
		EbfFileHT_close(fileht);

		if(0==fileht->ecode)
		{
			fileht->fp = fopen(filename, "rb+");
			if (fileht->fp != NULL)
			{
				EbfFileHT_close(fileht);
				EbfFileHT_create_htable(filename, 16);
				fileht->ecode=EbfTable_Put(filename,"/.ebf/info",0);
				fileht->ecode=EbfTable_Put(filename,"/.ebf/htable",ebfh.capacity_+ebfh.headersize);
			}
			else
			{
				fileht->ecode = 205;
			}
		}
	}
	else
	{
		printf("%s %s", "EBF error opening file:", filename);
		fileht->ecode = 205;
	}
	return fileht->ecode;
}


int EbfTable_InitSwap(const char* filename)
{
	int i;
	int64_t dims[8];
	EbfHeader ebfh;
	EbfFileHT fileht1;
	EbfFileHT* fileht;
	fileht=&fileht1;
	fileht->ecode=0;
	fileht->fp = fopen(filename, "wb");
	if (fileht->fp != NULL)
	{
		dims[0] = 5;
		EbfHeader_set(&ebfh, "/.ebf/info", 3, 8, "", "", 1, dims);
		for (i = 0; i < 5; ++i)
			fileht->hinfo[i] = 0;
		/* just making sure capacity_ is same as datasize, if in future _set does not do it*/
		ebfh.capacity_=EbfHeader_elements(&ebfh)*ebfh.datasize;
		ebfh.flag_swap=1;
		EbfHeader_write(&ebfh, fileht->fp);
		fileht->ecode=ebfh.ecode;
		ebfutils_SwapEndian(&(fileht->hinfo[0]),sizeof(fileht->hinfo[0]),EbfHeader_elements(&ebfh));
		if(fwrite(fileht->hinfo,ebfh.capacity_, 1, fileht->fp) != 1)
			fileht->ecode=206;
		EbfFileHT_close(fileht);

		if(0==fileht->ecode)
		{
			fileht->fp = fopen(filename, "rb+");
			if (fileht->fp != NULL)
			{
				EbfFileHT_close(fileht);
				EbfFileHT_create_htable(filename, 8);
				fileht->ecode=EbfTable_Put(filename,"/.ebf/info",0);
				fileht->ecode=EbfTable_Put(filename,"/.ebf/htable",ebfh.capacity_+ebfh.headersize);
			}
			else
			{
				fileht->ecode = 205;
			}
		}
	}
	else
	{
		printf("%s %s", "EBF error opening file:", filename);
		fileht->ecode = 205;
	}
	return fileht->ecode;
}


/*
void EbfTable_SwapEndian(const char* filename, int* ecode)
{
	EbfKeyVal* keyval=NULL;
	FILE* fp1;
	FILE* fout;
	EbfHeader header1,header2;
	char filename_out[1000];
	int flagswap,i,nsize;
	if((strlen(filename)+5+1)>1000)
	{
		*ecode=5;
	}
	else
	{
		strncpy(filename_out,filename);
		if(strncmp(filename_out[strlen(filename_out)-5],".ebf",4)==0)
		{
			filename_out[strlen(filename_out)-5]=0;
			strcat(filename_out,"_swap.ebf");
		}

	    fp1 = fopen(filename, 'rb');
	    fseek(fp1,0,SEEK_END);
	    filesize=ftell(fp1)
	    fseek(fp1,0,SEEK_SET);
	    EbfHeader_read(header1,fp1,ecode);
	    fseek(fp1,0,SEEK_SET);
	    if (header1.flag_swap == 0)
	    {
	        flagswap=1;
	        EbfTable_InitSwap(filename_out,ecode);
	    }
	    else
	    {
	        flagswap=0;
	        EbfTable_Init(filename_out);
	    }

	    fout = fopen(filename_out, 'rb+');
	    fseek(fout,0,SEEK_END);
	    while (fp1.Tell() < filesize)
	    {
	        EbfHeader_read(header1,fp1);
	        loc=ftell(fp1);
	        if (header1.datatype == 8)
	        {
	        	*ecode=11;
	        	printf("%s \n","Ebf error: structures cannot be swapped");
	        }
	        else
	        {
	            dth1 = TypeManager.itos_s(header1.datatype)
	        }
	        dblock_size=header1.elements()*header1.datasize

	        if ((strncmp(header1.name,'/.ebf/',6)!=0) && (strncmp(header1.name,'/.tr/',5)==False))
	        {
	            data = numpy.fromstring(fp1.Read(dblock_size), dtype = dth1);
	            if ((header1.flagswap == 1) && (flagswap==0))
	            {
	                data = data.byteswap(True);
	            }
	            if((header1.flagswap == 0) && (flagswap==1))
	            {
	                data = data.byteswap(True);
	            }
	            data = data.reshape(header1.getshape());
	            print header1.name,header1.getshape(), data.shape;
	            header2.create(header1.name,data,header1.dataunit,header1.sdef);
	            print header2.name,header2.getshape(), data.dtype;
	            if(flagswap == 1)
	            {
	                header2.flagswap = 1;
	            }

	            EbfHeader_write(header2,fout);
	            fwrite(data,dblock_size,1,fout)
	        }

	        fseek(fp1,loc+header1.capacity(), SEEK_SET);
	    }

	    filesize1=ftell(fp1);
	    fclose(fp1);
	    fclose(fout);
	    if(filesize1 != filesize)
	        print("%s \n","EBFCorrupt")
	    for(int i=0;i<nsize;++i)
	    {
	    	EbfTable_Put(filename_out,keyval[i].key, keyval[i].value,ecode);
	    }

	}

}
*/

int EbfTable_Remove(const char* filename, const char* keyin)
{
	int64_t loc,locp,hindex=0;
	EbfFileHTNode node,nodep;
	char key1[EBF_KEYNAME_MAXSIZE];
	char key[EBF_KEYNAME_MAXSIZE];
	size_t i;
	EbfFileHT fileht1;
	EbfFileHT* fileht;
	fileht=&fileht1;
	fileht->ecode=0;

	if(ebfstrcpy(key,EBF_KEYNAME_MAXSIZE,keyin)!=0)
		fileht->ecode=207;

	for (i = 0; i < strlen(key); ++i)
		key[i] = tolower(key[i]);


	if(fileht->ecode==0)
	{
		EbfFileHT_setup(fileht, filename,"rb+");
		hindex=EbfLtHash(key,fileht->header.htcapacity);
		EbfFileHT_read_hvalue(fileht, hindex,&loc);
	}

	if ((loc != 0)&&(fileht->ecode==0))
	{
		EbfFileHT_read_node(fileht, loc,&node);
		EbfFileHT_read_key(fileht, node.keyloc,node.keysize,key1);
		if(strcmp(key1,key)==0)
		{
			int64_t temp=0;
			if(node.next != -1)
				temp=node.next;
			EbfFileHT_write_hvalue(fileht, hindex,&temp);
		}
		else
		{
			locp=-1;
			while((node.next != -1)&&(strcmp(key1,key)!=0)&&(fileht->ecode==0))
			{
				locp=loc;
				nodep=node;
				loc=node.next;
				EbfFileHT_read_node(fileht, nodep.next,&node);
				EbfFileHT_read_key(fileht, node.keyloc,node.keysize,key1);
				if(node.next==nodep.next)
					fileht->ecode=217;
			}

			if(locp<0)
			{
				printf("%s \n","Ebf error: in remove(), something wrong with htable link list");
			}
			else
			{
				if(strcmp(key1,key)==0)
				{
					nodep.next=node.next;
					EbfFileHT_write_node(fileht, locp,&nodep);
				}
			}
		}
	}
	EbfFileHT_close(fileht);
	return fileht->ecode;

}


int Ebf_Rename(const char* filename, const char* oldkey, const char* newkey)
{
	EbfHeader ebfh;
	int64_t loc,loc1;
	int nexti=-1;
	size_t i;
	char mystr[EBF_KEYNAME_MAXSIZE+100];
	EbfFileHT fileht1;
	EbfFileHT *fileht;
	fileht=&fileht1;
	fileht->ecode=0;


	if((strlen(newkey)>=EBF_KEYNAME_MAXSIZE)||(strlen(oldkey)>=EBF_KEYNAME_MAXSIZE))
	{
		fileht->ecode=207;
		printf("%s \n","Ebf error Ebf_rename() strlen(newkey,oldkey)> max_size");
	}
	else
	{
		mystr[0]=0;
		if(strlen(newkey)==0)
		{
			do
			{
				nexti++;
				sprintf(mystr,"%s%s%s%d","/.tr",oldkey,".t",nexti);
				if(nexti>100000000)
				{
					fileht->ecode=220;
					break;
				}
			}
			while(EbfTable_Get(filename,mystr,&(fileht->ecode))!=-1);
		}
		else
		{
			if(ebfstrcpy(mystr,EBF_KEYNAME_MAXSIZE,newkey)!=0)
				fileht->ecode=207;
		}
	}

	if(strcmp(oldkey,mystr)==0)
	{
		fileht->ecode=219;
		printf("%s \n","Ebf error in Ebf_rename() oldkey == newkey");
	}


	if(fileht->ecode==0)
	{
		for(i=0;i<strlen(mystr);++i)
			mystr[i]=tolower(mystr[i]);

		loc=-1;
		if (fileht->ecode == 0)
			loc=EbfTable_Get(filename, oldkey,&(fileht->ecode));

		loc1=0;
		if (fileht->ecode == 0)
			loc1=EbfTable_Get(filename, mystr,&(fileht->ecode));

		EbfFileHT_setup(fileht, filename,"rb+");
		if(loc1>=0)
		{
			fileht->ecode=220;
			printf("%s \n","Ebf error: key with suggested name already exists");
		}
		if(loc<0)
		{
			printf("%s \n","Ebf error: key not found");
		}

		if((loc >= 0)&&(fileht->ecode==0)&&(loc1<0))
		{
			fseek(fileht->fp, loc, SEEK_SET);
			EbfHeader_read(&ebfh, fileht->fp);
			fseek(fileht->fp, loc, SEEK_SET);
			if (ebfh.ecode == 0)
				EbfHeader_rename(&ebfh, fileht->fp, mystr);
			EbfFileHT_close(fileht);
			if (ebfh.ecode == 0)
				ebfh.ecode=EbfTable_Remove(filename, oldkey);
			if (ebfh.ecode == 0)
				ebfh.ecode=EbfTable_Put(filename, mystr, loc);
			fileht->ecode = ebfh.ecode;
		}
		else
		{
			EbfFileHT_close(fileht);
		}
	}
	return fileht->ecode;
}



static void Ebf_GetKeyVal(const char* filename, EbfKeyVal **ebfh_list,
		size_t *hl_size, size_t* hl_capacity, int* ecode)
/* void Ebf_FileHT_getKeyValHT(const char* filename, char** key, int64_t** value,int *nsize1,int* ncapacity1) */
{
	FILE* fd = NULL;
	EbfHeader ebfh;
	EbfKeyVal *temp_p;
	size_t loc1;
	int64_t loc;
	size_t ncapacity;

	if ((*ebfh_list) != NULL)
	{
		printf(
				"%s \n",
				"EBF error: from Ebf_FileHT_getKeyVal() -> input headerlist should be NULL");
		*ecode = 8;
	}
	else
	{
		*hl_size=0;
		*hl_capacity=0;
		fd = fopen(filename, "rb");
		if (fd != NULL)
		{
			fseek(fd, 0, SEEK_END);
			loc = ftell(fd);
			fseek(fd, 0, SEEK_SET);
			while (ftell(fd) < loc)
			{
				loc1=ftell(fd);
				if((*hl_size)==(*hl_capacity))
				{
					ncapacity=(*hl_capacity)*2;
					if (ncapacity==0)
						ncapacity=16;
					temp_p = (EbfKeyVal*) realloc(*ebfh_list,ncapacity * sizeof(EbfKeyVal));
					if (temp_p != NULL)
					{
						*ebfh_list = temp_p;
						*hl_capacity=ncapacity;
					}
					else
					{
						*ecode=4;
						break;
					}
				}

				EbfHeader_read(&ebfh, fd);
				if (ebfh.ecode != 0)
				{
					*ecode = ebfh.ecode;
					break;
				}

				/* if size incorrect then not added to list. Hence no error reported or else it will terminate*/
				if(ebfstrcpy((*ebfh_list)[*hl_size].key,EBF_KEYNAME_MAXSIZE,ebfh.name)==0)
				{
					(*ebfh_list)[*hl_size].value=loc1;
					*hl_size=(*hl_size)+1;
				}
				else
				{
					printf("%s \n","Ebf Error: key size too large");
				}
				fseek(fd, ebfh.capacity_, SEEK_CUR);

			}
			fclose(fd);
		}
	}
}


static  void ebfutils_transferInSteps(FILE* fdin, FILE* fdout, int64_t size,int* ecode)
{
	size_t nread = 0, nwritten = 0;
	int64_t i = 0;
	int64_t step = 10000000;
	char* data;
	*ecode=0;
	data = (char *) malloc(step);
	if (data != NULL)
	{
		while ((i + step) <= size)
		{
			nread += fread(data, step, 1, fdin);
			nwritten += fwrite(data, step, 1, fdout);
			i += step;
		}
		step = size - i;
		if (step > 0)
		{
			nread += fread(data, step, 1, fdin);
			nwritten += fwrite(data, step, 1, fdout);
			i += step;
		}
		free(data);
		if (nread == nwritten)
			*ecode=0;
		else
			*ecode=6;
	}
	else
	{
		fclose(fdin);
		fclose(fdout);
		printf("%s \n","Ebf Error from transferInSteps()- Malloc Failure");
		*ecode=4;
	}

}

int64_t Ebf_GetCheckSum(const char* filename,int *ecode)
{
	EbfFileHT fileht1;
	EbfFileHT *fileht;
	fileht=&fileht1;
	fileht->ecode=0;
	EbfFileHT_setup(fileht, filename,"rb");
	EbfFileHT_read_hinfo(fileht);
	EbfFileHT_close(fileht);
	*ecode=fileht->ecode;
	return fileht->hinfo[0];
}






int Ebf_Copy(const char* infile, const char* outfile, const char* mode1,
		const char* dataname)
{
	size_t lh_size,lh_capacity;
	int ecode;
	char *c1=NULL;
	char *pch;
	EbfKeyVal* ebfh_list=NULL;
	FILE* fdin = NULL;
	FILE* fdout = NULL;
	char smode[10];
	int outexists=0;
	size_t i = 0;
	int64_t dataSize = 0;
	EbfHeader ebfh ;

	ecode=0;
	if (strcmp(infile, outfile) == 0)
	{
		printf("%s \n","Ebf Error from ebfcopy-Input and outout files should be different");
		ecode=1;
	}

	if ((strcmp(mode1, "w") != 0) && (strcmp(mode1, "a") != 0))
	{
		printf("%s \n","Ebf Error from ebfcopy():Mode must be w or a");
		ecode=1;
	}

	fdout = fopen(outfile, "rb");
	if (fdout != NULL)
	{
		fclose(fdout);
		fdout=NULL;
		outexists=1;
	}
	else
	{
		outexists=0;
	}



	/* no file exists is same as w mode */
	if ((strcmp(mode1, "w") == 0)||(outexists == 0))
		ecode=EbfTable_Init(outfile);

	/* for append r+ is needed to come back and overwrite the header if needed in ebf.Close() */
	strcpy(smode, "rb+");

	if (strlen(dataname)== 0)
	{
		lh_size=0;
		lh_capacity=0;
		Ebf_GetKeyVal(infile,&ebfh_list,&lh_size,&lh_capacity,&ecode);
	}
	else
	{
		c1=ebfstrdup(dataname);
		pch = strtok (c1,", ");
		lh_capacity=0;
		while (pch != NULL)
		{
			pch = strtok (NULL,", ");
			lh_capacity++;
		}
		lh_size=0;
		ebfh_list = (EbfKeyVal*) malloc(sizeof(EbfKeyVal) * lh_capacity);
		if (ebfh_list != NULL)
		{
			strcpy(c1,dataname);
			pch = strtok (c1,", ");
			while (pch != NULL)
			{
				if(ebfstrcpy(ebfh_list[lh_size].key,EBF_KEYNAME_MAXSIZE,pch)!=0)
				{
					ecode=7;
					break;
				}

				ebfh_list[lh_size].value=EbfTable_Get(infile,pch,&ecode);
				if(ebfh_list[lh_size].value<0)
				{
					printf("%s %s \n","dataname=",pch);
					printf("%s \n","Ebf Error from ebfcopy -data object not present in infile");
					ecode=9;
					break;
				}
				pch = strtok (NULL,", ");
				lh_size++;
			}
		}
		else
		{
			ecode = 4;
		}
		free(c1);
	}


	for(i=0;i<lh_size;++i)
	{
		if ((strncmp(ebfh_list[i].key,"/.ebf/",6) != 0)&&(strncmp(ebfh_list[i].key,"/.tr/",5) != 0))
		{
			if(EbfTable_Get(outfile,ebfh_list[i].key,&ecode)>0)
			{
				printf("%s %s \n","dataname=",ebfh_list[i].key);
				printf("%s \n","Ebf Error from ebfcopy -data object already present in outfile");
				ecode=9;
				break;
			}
		}
	}



	/* Open files */
	fdin = fopen(infile, "rb");
	if (fdin == NULL)
	{
		printf("%s %s \n","Ebf Error from ebfcopy: cannot open out file" , infile);
		ecode=1;
	}


	fdout = fopen(outfile, smode);
	if (fdout == NULL)
	{
		printf("%s %s \n","Ebf Error from ebfcopy: cannot open out file" , infile);
		ecode=1;
		printf("%s \n", outfile);
	}

	/* compute datasize and position files for I/O */
	if ((ecode) == 0)
	{
		fseek(fdout, 0, SEEK_END);
		for(i=0;i<lh_size;++i)
		{
			fseek(fdin, ebfh_list[i].value, SEEK_SET);
			ebfh_list[i].value=-1;
			EbfHeader_read(&ebfh, fdin);
			if ((0==(ebfh.ecode))&&(strncmp(ebfh.name,"/.ebf/",6) != 0)&&(strncmp(ebfh.name,"/.tr/",5) != 0))
			{
				ebfh_list[i].value=ftell(fdout);
				ebfh.version[0]=1;
				ebfh.version[1]=1;
				ebfh.capacity_=EbfHeader_elements(&ebfh)*ebfh.datasize;
				dataSize = ebfh.capacity_;
				EbfHeader_write(&ebfh, fdout);
				if(ebfh.ecode==0)
				{
					ebfutils_transferInSteps(fdin, fdout, dataSize,&ecode);
				}
				else
				{
					printf("%s \n","Error while writing");
					lh_size=i+1;
				}
			}
		}
	}
	else
	{
		lh_size=0;
	}


if (fdout!=NULL)
	fclose(fdout);
if (fdin!=NULL)
	fclose(fdin);

for(i=0;i<lh_size;++i)
{
	if(ebfh_list[i].value>=0)
		ecode=EbfTable_Put(outfile,ebfh_list[i].key,ebfh_list[i].value);
}

if(ebfh_list!=NULL)
	free(ebfh_list);



#ifdef EBFLANGOPTCPP
if((ecode)!=0)
{
/*	printf("%s \n","EBF Error from ebfcopy()"); */
	throw std::logic_error("EBF Error from ebfcopy()");
}
#endif

/*
	checksum=Ebf_GetCheckSum(outfile);
	sprintf(mystr,"%s%s%s %"PRId64"%s","(",dataname,",",checksum,") ");
	if(strlen(mystr) >= (EBF_KEYNAME_MAXSIZE+100))
		ebfutils_error("Ebf Error from EbfFile.Close()- mystr max size is too small");
	checksum=EbfCkHash(mystr,checksum);
	Ebf_updateCheckSum(efile->filename,checksum);

if ((strlen(dataname) == 0)&&((*ecode)==0))
{
	printf("%s %"PRId64" \n",infile, Ebf_GetCheckSum(infile,ecode));
	printf("%s %"PRId64" \n",outfile, Ebf_GetCheckSum(outfile,ecode));
}
*/

return ecode;
}


int Ebf_ContainsKey(const char* filename1, const char* dataname1,int* ecode)
{
	*ecode=0;
	if(EbfTable_Get(filename1,dataname1,ecode)>=0)
		return 1;
	else
		return 0;
}


int Ebf_ContainsKey_Info(const char * filename1, const char* dataname1,EbfDataInfo* dinfo)
{
	EbfHeader ebfh;
	int i;
	dinfo->ecode = 0;
	dinfo->datapos = -1;
	dinfo->rank = 0;
	dinfo->datasize = 0;
	dinfo->elements = 0;
	dinfo->headerpos = -1;
	for (i = 0; i < EBF_RANK_MAXSIZE; ++i)
		dinfo->dim[i] = 0;
	ebfh.ecode=0;
	dinfo->headerpos=EbfTable_Get_WithHeader(filename1,dataname1,&ebfh);
	dinfo->ecode = ebfh.ecode;

	if (dinfo->headerpos >= 0)
	{
		dinfo->datasize = ebfh.datasize;
		dinfo->datatype = ebfh.datatype;
		dinfo->rank = ebfh.rank;
		dinfo->datapos = ebfh.datapos;
		if (ebfh.rank > EBF_RANK_MAXSIZE)
		{
			dinfo->ecode = 11;
			printf("%s \n","Rank greater than 8 not allowed");
		}
		else
		{
			for (i = 0; i < ebfh.rank; ++i)
				dinfo->dim[i] = ebfh.dim[i];
			dinfo->elements = EbfHeader_elements(&ebfh);
			dinfo->datapos = ebfh.datapos;
			strncpy(dinfo->unit, ebfh.unit, EBF_KEYNAME_MAXSIZE);
			dinfo->unit[EBF_KEYNAME_MAXSIZE-1]=0;
		}
	}

	if((dinfo->ecode==0)&&(dinfo->headerpos>=0))
		return 1;
	else
		return 0;
}



#ifdef EBFLANGOPTCPP
void EbfTable_Print(const char* filename)
{
	EbfKeyVal* keyval=NULL;
	int i,nsize=0;
	int ncapacity=0;
	int ecode=0;
	EbfFileHT fileht1;
	EbfFileHT *fileht;
	fileht=&fileht1;
	std::cout<<std::setw(24)<<"filename:"<<filename<<std::endl;
	EbfFileHT_setup(fileht,filename,"rb");
	std::cout<<"\t"<<std::setw(24)<<"ecode="<<fileht->ecode<<std::endl;
	std::cout<<"\t"<<std::setw(24)<<"hecode="<<fileht->hecode<<std::endl;
	std::cout<<"\t"<<std::setw(24)<<"hinfo flagswap="<<fileht->ebfh_hinfo.flag_swap<<std::endl;
	std::cout<<"\t"<<std::setw(24)<<"flagswap="<<fileht->flagswap<<std::endl;
	std::cout<<"\t"<<std::setw(24)<<"hinfo="<<fileht->hinfo[0]<<" "<<fileht->hinfo[1]<<" "<<fileht->hinfo[2]<<" "<<fileht->hinfo[3]<<std::endl;
	if(fileht->ecode==0)
	{
		std::cout<<"header:"<<std::endl;;
		std::cout<<"\t"<<std::setw(24)<<"count="<<fileht->header.count<<std::endl;
		std::cout<<"\t"<<std::setw(24)<<"htcapapcity="<<fileht->header.htcapacity<<std::endl;
		std::cout<<"\t"<<std::setw(24)<<"nodecapapcity="<<fileht->header.nodecapacity<<std::endl;
		std::cout<<"\t"<<std::setw(24)<<"keycapapcity="<<fileht->header.keycapacity<<std::endl;
		std::cout<<"\t"<<std::setw(24)<<"keyposcur="<<fileht->header.keyposcur<<std::endl;
	}
	EbfFileHT_close(fileht);

	if(fileht->ecode==0)
	{
		std::cout<<"Key Value List:"<<std::endl;
		EbfTable_getKeyValHT(filename,&keyval,&nsize,&ncapacity,&ecode);
		if(ecode==0)
		{
			std::cout<<"nsize="<<std::setw(24)<<nsize<<std::endl;
			for(i=0;i<nsize;++i)
				std::cout<<std::setw(24)<<keyval[i].key<<" -> "<<keyval[i].value<<std::endl;
		}
		else
		{
			std::cout<<"Ebf Error, ecode!=0"<<std::endl;
		}
		if(keyval!=NULL)
		{
			free(keyval);
			keyval=NULL;
		}
	}
	else
	{
			std::cout<<"Ebf Error, ecode!=0"<<std::endl;
	}

}





}
}
#else
#ifdef EBFLANGOPTIDL
void EbfTable_Print(const char* filename)
{
}
#else
void EbfTable_Print(const char* filename)
{
	EbfKeyVal* keyval=NULL;
	int i,nsize=0;
	int ncapacity=0;
	int ecode=0;
	EbfFileHT fileht1;
	EbfFileHT *fileht;
	fileht=&fileht1;
	printf("%-24s %s \n","filename:",filename);
	EbfFileHT_setup(fileht,filename,"rb");
	printf("\t %-24s %d \n","ecode=",fileht->ecode);
	printf("\t %-24s %d \n","hecode=",fileht->hecode);
	printf("\t %-24s %d \n","hinfo_flagswap=",fileht->ebfh_hinfo.flag_swap);
	printf("\t %-24s %d \n","flagswap=",fileht->flagswap);
	printf("\t %-24s %"PRId64" %"PRId64" %"PRId64" \n","hinfo=",fileht->hinfo[0],fileht->hinfo[1],fileht->hinfo[2]);
	if(fileht->ecode==0)
	{
		printf("%s \n","header:");
		printf("\t %-24s %"PRId64"\n","count=",fileht->header.count);
		printf("\t %-24s %"PRId64"\n","htcapapcity=",fileht->header.htcapacity);
		printf("\t %-24s %"PRId64"\n","nodecapapcity=",fileht->header.nodecapacity);
		printf("\t %-24s %"PRId64"\n","keycapapcity=",fileht->header.keycapacity);
		printf("\t %-24s %"PRId64"\n","keyposcur=",fileht->header.keyposcur);
	}
	EbfFileHT_close(fileht);

	if(fileht->ecode==0)
	{
		printf("%s \n","Key Value List:");
		EbfTable_getKeyValHT(filename,&keyval,&nsize,&ncapacity,&ecode);
		if(ecode==0)
		{
			printf("%s %24d \n","nsize=",nsize);
			for(i=0;i<nsize;++i)
				printf("%s %"PRId64"\n",keyval[i].key,keyval[i].value);
		}
		else
		{
			printf("%s \n","Ebf Error, ecode!=0");
		}
		if(keyval!=NULL)
		{
			free(keyval);
			keyval=NULL;
		}
	}
	else
	{
		printf("%s \n","Ebf Error, ecode!=0");
	}

}
#endif
#endif


/*-------------------------------------------------------------*/


