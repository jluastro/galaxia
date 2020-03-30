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

#ifndef COT_H_
#define COT_H_
#include "Matrix.h"

class Cot
{
public:
	Cot();
	static void print() {TM.print();}
	static void xyz_to_lbr(double *PosC,int dim=3);
	static void xyz_to_lbr(double *PosC,double *PosS,int dim=3);
	static void lbr_to_xyz(double *PosS,int dim=3);
	static void lbr_to_xyz(double *PosS,double *PosC,int dim=3);
	static void lb_to_radec(double l,double b,double &ra,double& dec);

	static void xyz_to_lzr(double *PosC,int dim=3);
	static void xyz_to_lzr(double *PosC,double *PosS,int dim=3);
	static void lzr_to_xyz(double *PosS,int dim=3);
	static void lzr_to_xyz(double *PosS,double *PosC,int dim=3);
	static void rotate(double *PosC);

	static void setTM(double l1,double b1);
	static void setTM_lzr(double l1,double b1);
	static const Matrix<double>& RotationMatrix(double* Pos,double theta);
	static const Matrix<double>& EulerMatrix(double phi,double theta, double psi);
	static const Matrix<double>& YPRMatrix(double phi,double theta, double psi);
	static const Matrix<double>& RotationMatrix(int axis,double theta);
private:
	static double l,b,r,sin_b,cos_b,sin_l,cos_l;
	static std::vector<double> x;
	static Matrix<double> TM;
	static Matrix<double> TR;
	~Cot();
};

#endif /* COT_H_ */
