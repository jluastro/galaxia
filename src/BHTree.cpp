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

#include<algorithm>
#include "BHTree.h"

//static double sigma=20.0;
//static double fd=sigma*sqrt(2*3.14159);
//

//double function1(double* Pos)
//{
//	double rad2=radius2(Pos);
////	cout<<Pos[0]<<" "<<rad<<endl;
//	return 1e8*exp(-(rad2/(2*sigma*sigma)))/(fd*fd*fd);
//}


//double function1(double* Pos)
//{
//	double rho0,d0,hr1,hr2,a2,eps;
//	hr1=2.530; hr2=1.320; eps=0.0268,d0=0.032,rho0=7.9e6;
//	a2=radius2(Pos)+(Pos[2]*Pos[2])/eps;
////	cout<<Pos[0]<<" "<<rad<<endl;
//	return (rho0/d0)*(exp(-sqrt(0.25+a2/(hr1*hr1)))-exp(-sqrt(0.25+a2/(hr2*hr2))));
//}

//double function1(double* Pos)
//{
//	double rho0,d0,hr1,hr2,a2,eps;
//	hr1=3.5; hr2=0.2; eps=0.0268,d0=0.0881,rho0=7.9e6;
//	a2=sqrt(Pos[0]*Pos[0]+Pos[1]*Pos[1]);
//	return (rho0/d0)*exp(-(a2/hr1+fabs(Pos[2])/hr2));
//}



//BHNode::BHNode()
//{
////	left=NULL;
//}

BHNode::~BHNode()
{
}

BHTree::BHTree(Population &Pop,const double* posC,int warpFlareOn1,int option,string inputDir)
{
	initialize(Pop,posC,warpFlareOn1,option,inputDir);
}

void BHTree::initialize(Population &Pop, const double* posC,int warpFlareOn1,int option,string inputDir)
{
	string fname;

//	sprintf(fname, "%s%s%i%s","data_bhtree/", outputFile, Pop.ID, ".ebf");

	if((warpFlareOn1)>0)
	{
		stringstream sout;
		if(Pop.optionE==1)
			sout<<inputDir<<"BHTree-2.2/bhtree_with_wf/"<<outputFile<<Pop.ID<<"_E1.ebf";
		else
			sout<<inputDir<<"BHTree-2.2/bhtree_with_wf/"<<outputFile<<Pop.ID<<"_E0.ebf";
		fname=sout.str();
//		sprintf(fname, "%s%s%s%i%s",inputDir,"/bhtree_with_wf/", outputFile, Pop.ID, ".ebf");
	}
	else
	{
		stringstream sout;
		if(Pop.optionE==1)
		sout<<inputDir<<"BHTree-2.2/bhtree_no_wf/"<<outputFile<<Pop.ID<<"_E1.ebf";
		else
		sout<<inputDir<<"BHTree-2.2/bhtree_no_wf/"<<outputFile<<Pop.ID<<"_E0.ebf";
		fname=sout.str();
//		sprintf(fname, "%s%s%s%i%s",inputDir,"/bhtree_no_wf/", outputFile, Pop.ID, ".ebf");
	}
//	FILE* fd;
//	if ((fd = fopen(fname, "r")))
//	{
//		load(fname);
//	}
//	else

//	if (Pop.ID==8)
	if (option==0)
	{
		cout << "Making Tree" << endl;
		rootv.resize(10000000);
		root = &rootv[0];
		{
			vector<double> x(3, 0.0);
			root[0].setUpRoot(&x[0], 300.0, Pop);
		}
		root[0].computeMassRefined(1.0, 3, Pop);
		//	cout<<root[0].mass<<endl;

		vector<BHNode*> nodeStack;
		nodeStack.push_back(&root[0]);
		BHNode* lNode = &root[0];
		BHNode* cNode;
		while (nodeStack.size() > 0)
		{
			cNode = nodeStack.back();
			nodeStack.pop_back();
			if (cNode->splittable(Pop.mass_split,Pop.l_split))
			{
//				cout<<(lNode-root)<<" "<<rootv.size()<<endl;
				assert(size_t(lNode-root+8)<rootv.size());
				cNode->split(lNode, Pop);
				for (int i = 0; i < 8; ++i)
				{
					lNode++;
					nodeStack.push_back(lNode);
				}
			}
			else
				leafs.push_back(cNode);
			//		cout << nodeStack.size() << " " << " " << lNode - root<<" "<<cNode->mass<<" "<<cNode->delta_mass/cNode->mass<<" "<<cNode->l<< endl;
			//		exit(1);
		}
		rootv.resize((lNode - root) + 1);
		assert(root==&rootv[0]);

		cout<<" Total Nodes="<<rootv.size()<<" leafs="<<leafs.size()<<endl;
		cout<<"Completed %  <";
		double temp = 0.0;
		double temp1 = 0.0;
		double temp2 = 0.0;
		for (size_t i = 0; i < leafs.size(); ++i)
		{
			leafs[i]->computeMassRefined(1.0, 5, Pop);
			temp += leafs[i]->mass;
			if(leafs[i]->mass>50.0)
			{
				temp1+=leafs[i]->rho_max*leafs[i]->mass/leafs[i]->rho_min;
				temp2+=leafs[i]->mass;
			}
			if((i%(leafs.size()/10))==0)
				cout<<(i*10*10)/leafs.size()+1<<".."; cout<<flush;
		}
		cout<<">"<<endl;
		cout <<"Mass=" << temp << " Leafs=" << leafs.size()<<" Total Nodes="<<rootv.size()<<" "<<temp2<<" "<<temp1/temp2<< endl;

		// check
//		{
//			double temp = 0.0;
//			for (size_t i = 0; i < leafs.size(); i += 1000)
//			{
//				temp += leafs[i]->mass;
//				leafs[i]->calculateDmod(posC);
//				cout << i << " " << leafs[i]->mass << " " << temp << endl;
//			}
//		}
		save(fname);
	}
	if (option==1)
	{
		FILE* fd;
		if ((fd = fopen(fname.c_str(), "r")))
		{
			fclose(fd);
			load(fname);

// check
//			{
//				double temp = 0.0;
//				for (size_t i = 0; i < leafs.size(); i += 1000)
//				{
//					temp += leafs[i]->mass;
//					leafs[i]->calculateDmod(posC);
//					cout << i << " " << leafs[i]->mass << " " << temp << endl;
//				}
//			}
//			exit(1);
			double temp = 0.0;
			for (size_t i = 0; i < leafs.size(); ++i)
			{
				temp += leafs[i]->mass;
				leafs[i]->calculateDmod(posC);
//				cout<<i<<" "<<leafs[i]->mass<<" "<<temp<<endl;
			 }
//			cout<<left<<setw(36)<<"Mass="<<setw(24)<<temp<<endl;
//			cout<<left<<setw(36)<<"Leafs="<<setw(24)<<leafs.size()<<endl;
			sort(&leafs[0], &leafs[leafs.size()], nodep_less_than());
		}
		else
		{
			cout<<"BHtree file not found."<<fname<<endl;
			cout<<"Please run the setup with -s option. "<<endl;
			cout<<"For info check with --help option"<<endl;
			exit(1);
		}
	}
}

BHTree::~BHTree()
{
	// TODO Auto-generated destructor stub
}


//		if(leafs[i]->mass>1e2)
//		{
//			tempv1.push_back(leafs[i]->mass);
//			tempv2.push_back(leafs[i]->l);
//		}
//	}
//	cout<<"Mass="<<temp<<" "<<temp1/temp<<"Leafs="<<leafs.size()<<endl;
//	Stats(tempv1.begin(),tempv1.end());
//	Stats(tempv2.begin(),tempv2.end());



//	BHNode* root1=root;
//	vector<BHNode> rootv1=rootv;
//	vector<BHNode*> leafs1=leafs;

//	load(fname);

//	for(size_t i=0;i<rootv1.size();++i)
//		if(rootv1[i].left!=NULL)
//			rootv1[i].left=&rootv[0]+(rootv1[i].left-root1);
//	for(size_t i=0;i<rootv1.size();++i)
//		assert(rootv[i]==rootv1[i]);
//	for(size_t i=0;i<leafs1.size();++i)
//		assert((leafs[i]-root)==(leafs1[i]-root1));
