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

#ifndef EBFVECTOR_HPP_
#define EBFVECTOR_HPP_

#include "ebf.hpp"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>


namespace ebf
{

template<class T>
class EbfVector
{
public:
	EbfVector(){writeOn=0;};
	EbfVector(const std::string &fname, const std::string& dataname, std::string mode="r")
	{
		writeOn=0;
		init(fname,dataname,mode);
	}
	void init(const std::string &fname,const std::string& dataname,std::string mode="r")
	{
		ecode=0;
		if(writeOn!=0)
		{
//			std::cout<<"EbfVector: was Open in Write mode, hence cannot be reinitialized"<<std::endl;
			close();
			throw std::logic_error("Ebf Error from ebfvector-was open in write mode, hence cannot be reinitialized");
		}
		debug = 0;
		bufsize = 1000;
		fname_ = fname;
		dataname_ = dataname;
		datatype_ = TypeT<T>();
		if (mode=="w")
			writeOn = 1;
		else if (mode=="r")
			writeOn = 0;
		else
		{
			throw std::logic_error("Ebf Error from ebfvector-Must be read or write");
		}
		setup();
	}
	~EbfVector()
	{
		close();
	}
	inline T& operator [](int64_t i)
	{
		if ((i < offset) || (i >= (int64_t(value.size()) + offset)))
			loadNext(i);
		return value[i - offset];
	}
	inline T& operator ()(int64_t i1,int64_t i2)
	{
		int i=efile.Index(i1,i2);
		if ((i < offset) || (i >= (int64_t(value.size()) + offset)))
			loadNext(i);
		return value[i - offset];
	}
	inline T& operator ()(int64_t i1,int64_t i2,int64_t i3)
	{
		int i=efile.Index(i1,i2,i3);
		if ((i < offset) || (i >= (int64_t(value.size()) + offset)))
			loadNext(i);
		return value[i - offset];
	}
	void push_back(const T& a)
	{
		if (writeOn == 2)
		{
			value.push_back(a);
			size_++;
			if (int64_t(value.size()) == bufsize)
				tempflush();
		}
		else
		{
//			std::cout
//					<< "ERROR: in EbfVector, Write has been completed, or it is Read only"
//					<< std::endl;
			close();
			throw std::logic_error("Ebf Error from ebfvector:write has been completed, or it is read only");
//			exit(1);
		}
	}
	void flush()
	{
		if (writeOn == 2)
		{
			tempflush();
			efile.Close();
			writeOn = 3;
		}
	}
	size_t size()
	{
		return size_;
	}
	int rank()
	{
		return rank_;
	}
	const int64_t& dim(int i)
	{
		return dims_[i];
	}
	const std::string& getFileName()
	{
		return fname_;
	}
	void saveTo(const std::string& filename,const std::string& mode)
	{
		if(writeOn==2)
			flush();

		if(writeOn==3)
		{
			if(filename==fname_)
			{
//				std::cout<<"output file should not be same as input file"<<std::endl;
				close();
				throw std::logic_error("Ebf Error from ebfvector-output file should not be same as input file");
			}
			ecode=ebfc::Ebf_Copy(fname_.c_str(),filename.c_str(),mode.c_str(),"");
			close();
			writeOn=4;
		}
		else
		{
//			std::cout<<"the data is not flushed or it is in Read mode"<<std::endl;
			close();
			throw std::logic_error("Ebf Error from ebfvector-the data is not flushed or it is in read mode");
		}
	}
	const std::string& getDataName()
	{
		return dataname_;
	}
	void printStats()
	{
		double xmin,xmax,xmean,xstd;
		double temp=this->operator [](0);
		xmin=temp;
		xmax=temp;
		xmean=temp;
		xstd=temp*temp;
		for(int i=1;i<size_;++i)
		{
			temp=this->operator [](i);
			if (xmin > temp)
				xmin = temp;
			if (xmax < temp)
				xmax = temp;
			xmean += temp;
			xstd += temp*temp;
		}
		xmean=xmean/size_;
		xstd=xstd/size_;
		std::cout<<std::left<<std::setw(14)<<dataname_;
		std::cout<<std::setw(12)<<xmin<<std::setw(12)<<xmax<<std::setw(12)<<xmean<<std::setw(12)<<sqrt(xstd-xmean*xmean);
		std::cout << "(";
		for (int i = 0; i < efile.rank() - 1; ++i)
			std::cout << efile.dim(i) << ",";
		std::cout << efile.dim(efile.rank() - 1) << ")" << std::endl;


	}
//	EbfVector( const EbfVector& other)
//	{
//		debug=other.debug;
//		writeOn=other.writeOn;
//		rank=other.rank_;
//		datatype=other.datatype_;
//		throw std::runtime_error("Copying an EbfVector not allowed");
//	}
private:
	// un-muatble
	bool debug;
	int writeOn,ecode;
	int rank_, datatype_;
	std::vector<int64_t> dims_;
	std::string dataname_, fname_;
	int64_t bufsize;
	//----- mutable
	std::vector<T> value;
	int64_t offset, size_;
	EbfFile efile;
	EbfVector& operator=( const EbfVector& rhs )
	{
		throw std::runtime_error("assignment operator on an ebfvector not allowed");
	}
	void tempflush()
	{
		if (value.size() > 0)
		{
			if (debug)
				std::cout << "Writing flush" << std::endl;
			offset += value.size();
			efile.Write(&value[0], value.size());
			value.clear();
		}
	}
	void close()
	{
		efile.Close();
		if ((writeOn == 2) || (writeOn == 3))
		{
			if (ebfutils::FileExists(fname_) == 1)
				remove(fname_.c_str());
		}
	}
	void setup()
	{
		efile.Close();
		if (value.size() > 0)
			value.clear();
		offset = 0;

		if (writeOn == 1)
		{
			size_ = 0;
			if (ebfutils::FileExists(fname_) == 1)
			{
//				std::cout << "WARNING: File-" << fname_
//						<< " already in use or exists" << std::endl;
//				std::cout << "To avoid risk- delete the file and run again" << std::endl;
//				std::cout
//						<< "Also check if multiple instances of program using same file"
//						<< std::endl;
//				char c;
//				std::cout
//						<< "Do you understand the risk and want to overwrite? (Y/N):";
//				cin >> c;
//				if (c != 'Y')
//					exit(1);
			}
			efile.Open(fname_, dataname_, "w", TypeT<T>());
//			if (efile.getStatus() == 1)
				writeOn = 2;
		}
		else if (writeOn == 0)
		{
			efile.Open(fname_, dataname_);
			rank_ = efile.rank();
			size_ = efile.dim(0);
			dims_.push_back(efile.dim(0));
			for (int i = 1; i < efile.rank(); ++i)
			{
				if (efile.dim(i) > 0)
				{
					dims_.push_back(efile.dim(i));
					size_ *= efile.dim(i);
				}
			}
			loadNext(0);
//			if(efile.rank()==1)
//				dims_.push_back(1);
		}
		else
		{
//			std::cout << "Init not allowed with writeOn=" << writeOn << std::endl;
			close();
			throw std::logic_error("Ebf Error from ebfvector:Init not allowed with write on");
		}

//		if (efile.getStatus() != 1)
//		{
//			std::cout << "EbfVector: file and/or record in it not found-"<<std::endl;
//			std::cout<<"   "<< fname_ <<"["<<dataname_<<"]"<< std::endl;
//			Close();
//			throw std::runtime_error("Ebf Error from EbfVector:");
//		}
	}
	void loadNext(int64_t j)
	{
		if (writeOn != 0)
		{
//			std::cout << "EbfVector is not in Read mode" << std::endl;
			close();
			throw std::logic_error("Ebf Error from ebfvector-not in read mode");
		}
		if((j >= size_) || (j < 0))
		{
//				std::cout<<"ebfvevtor: out of bound Index for operator [] i="<<j<<" size="<<size_<<std::endl;
				throw std::range_error("Ebf Error from ebfvector:out of bound index for operator []");
		}
		if((offset > size_) && (offset < 0))
		{
//			std::cout<<"ebfvevtor: out of bound Index of var offset"<<std::endl;
			throw std::range_error("Ebf Error from ebfvector- out of bound index of var offset");
		}
		if (int64_t(j) < offset)
		{
			value.clear();
			offset = 0;
			efile.Seek(0);
//			setup();
		}

		while (j >= (int64_t(value.size()) + offset))
		{
			offset += value.size();
			if(j >= (int64_t(value.size()) + offset))
			{
				offset=j;
			}

			if ((offset + bufsize) > size_)
				value.resize(size_ - offset);
			else
				value.resize(bufsize);

			//-------
			if(TypeT<T>()==7)
			{
				efile.Seek(offset);
				efile.ReadAsStringV(value);
			}
			else
			{
				efile.Seek(offset);
				efile.Read(&value[0],value.size());
			}


			//--------
		}
	}

};


}


#endif /* EBFVECTOR_HPP_ */
