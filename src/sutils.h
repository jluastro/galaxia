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

#ifndef SUTILS_H_
#define SUTILS_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <ctime>

#include <stdint.h>
#include <inttypes.h>


using namespace std;

const double PI=3.14159265358979323846;

inline void certify(int x,const string& s="")
{
	if(x==0) {cout<<"CHECK FAILED: "<<s<<endl; exit(1);}
}
void stringSplit(const string &s,const char* delimiters,vector<string> &sv);

#endif /* SUTILS_H_ */
