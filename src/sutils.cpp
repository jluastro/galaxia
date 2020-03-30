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

#include "sutils.h"



bool fileExists(const char* filename)
{
	ifstream infile;
	infile.open (filename);
	if (infile.is_open())
	{
		infile.close();
		return 1;
	}
	else
	{
		return 0;
	}
}




void stringSplit(const char *c,const char* delimiters,vector<string> &sv)
{
	sv.clear();
//	vector<char> cv; 	cv.resize(strlen(c)+1);		char *c=&cv[0];
	char* c1= new char[strlen(c)+1];
	strcpy(c1,c);

	char * pch;
	pch = strtok (c1,delimiters);
	while (pch != NULL)
	{
		sv.push_back(pch);
	    pch = strtok (NULL,delimiters);
	}

	delete[] c1;
}


void stringSplit(const string &s,const char* delimiters,vector<string> &sv)
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


void stringSplit(stringstream &sout,const char* delimiters,vector<string> &sv)
{
	sv.clear();
	sout.seekg (0, ios::end);
	int length = sout.tellg();
	char* c1= new char[length+1];
	sout.seekg (0, ios::beg);
	sout.read (c1,length);
	sout.seekg (0, ios::beg);
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

//void stringstreamTovchar(stringstream &sout,vector<char> &cv)
//{
//	int64_t pos = sout.tellg();
//	sout.seekg (0, ios::end);
//	int length = sout.tellg();
//	cv.resize(length+1);
//	sout.seekg (0, ios::beg);
//	sout.read (&cv[0],length);
//	cv[length]=0;
//	sout.seekg (pos, ios::beg);
//}



