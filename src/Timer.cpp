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

#include "Timer.h"

Timer::Timer()
{
	start();
}

Timer::~Timer()
{
}

void Timer::start()
{
	startTime=clock();
	stopTime=clock();
}

void Timer::stop()
{
	stopTime=clock();
}


void Timer::resume()
{
	startTime=clock()-(stopTime-startTime);
	stopTime=clock();
}


double Timer::elapsedTime()
{
	return double(stopTime-startTime)/CLOCKS_PER_SEC;
}

double Timer::currentTime()
{
	return double(clock()-startTime)/CLOCKS_PER_SEC;
}


void Timer::print()
{
	stop();
	std::cout<<left<<setw(12)<<" Time "<<setw(12)<<elapsedTime()<<std::endl;
}
void Timer::print(const char *s)
{
	stop();
	std::cout<<left<<setw(36)<<s<<setw(12)<<elapsedTime()<<std::endl;
}
void Timer::print(const char *s,double n)
{
	stop();
	std::cout<<left<<setw(36)<<s<<setw(12)<<elapsedTime()<<"     (  "<<n*1e-6/elapsedTime()<<" Mops )"<<std::endl;
}

void Timer::printC()
{
	std::cout<<left<<setw(12)<<" Time "<<setw(12)<<currentTime()<<std::endl;
}
void Timer::printC(const char *s)
{
	std::cout<<left<<setw(36)<<s<<setw(12)<<currentTime()<<std::endl;
}
void Timer::printC(const char *s,double n)
{
	std::cout<<left<<setw(36)<<s<<setw(12)<<currentTime()<<"     (  "<<n*1e-6/currentTime()<<" Mops )"<<std::endl;
}

