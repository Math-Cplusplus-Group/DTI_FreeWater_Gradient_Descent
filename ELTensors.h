/** \file  ELTensors.h
\brief C++ source file initializing tensors.
Copyright 2016 by Andrew Colinet, Tomas Kojar

Redistribution and use in source and binary forms, with or without modification, are permitted provided
that the following conditions are met:

-# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

-# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

-# The name of the copyright holder may not be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELTENSORS_H
#define ELTENSORS_H
#include "ELInitialization.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>


class ELTensors {

private:
	//4-tensor of Cell X
	//for example CellX[0][0][1]=X(x,y,z+1) and CellX[0][0][2]=X(x,y,z-1)	
	std::vector<std::vector<std::vector<std::vector<double> > > > CellX;

	//mesh lengths
	double dx;
	double dy;
	double dz;
	//volume fraction
	double Volfn;
	ELInitialization Eli;
	std::vector<double> Ahat;

public:
	std::ofstream& myfile;

	/// Constructor.
	ELTensors(std::vector<std::vector<std::vector<std::vector<double> > > > xCellX,
		double xdx,
		double xdy,
		double xdz,
		double xVolfn,
		ELInitialization xEli,
		std::vector<double> xAhat,
		std::ofstream& xmyfile
		) : CellX(xCellX), dx(xdx), dy(xdy), dz(xdz), Volfn(xVolfn), Eli(xEli), Ahat(xAhat),myfile(xmyfile) { };


	/// Copy constructor.



	////Delete constructor.
	//~ELTensors();

												////Global variables needed for the EL equation

	double term1 = 0;
	double bval = Eli.bval;
	double  Awater = Eli.Awater;
	double alpha = Eli.alpha;
	int nuframesx = Eli.nframesx;
	int nuframesy = Eli.nframesy;
	int nuframesz = Eli.nframesz;
	int Graddirections = Eli.GradDirections;


													////Global tensors containers needed for the EL equation


	////They will be populated by the member functions



													////Tensors
	//// Iwasawa coordinates
	double w1 = CellX[0][0][0][3];
	double w2 = CellX[0][0][0][4];
	double w3 = CellX[0][0][0][5];
	double w4 = CellX[0][0][0][6];
	double 	w5 = CellX[0][0][0][7];
	double w6 = CellX[0][0][0][8];
	
	///Diffusion tensor	
	double Diff[3][3] = { { w1, w1*w4, w1*w5 }  ,
	{ w1*w4, w2 + w1*w4*w4,w1*w4*w5 + w2*w6 },
	{ w1*w5, w1*w4*w5 + w2*w6, w3 + w1*w5*w5 + w2*w6*w6 } };


	////Image metric h_ij
	std::vector<double> hmatrix ;
	
	////Induced Pullback metric gamma
	std::vector<std::vector<double>> IMmatrix;


	////Inverse of Induced Pullback metric gamma
	std::vector<std::vector<double>> InvIMmatrix;
	
	
	///The exponent (q_k^T) * Diffusion matrix * (q_k)
	std::vector<double> qDqsum;

	////Christoffel symbols
	std::vector<std::vector<std::vector<double> > > ChrisSym;




	
											///Derivatives of tensors

	////Derivative of embedding map X
	std::vector<std::vector<double>> DerivX;

	////Double Derivative of embedding map X
	std::vector<std::vector<std::vector<double> > > DDerivX;


	///Derivatives of Diffusion tensor
	double Dx1[3][3] = { { 1, w4, w5 }  ,
	{ w4, w4*w4, w4*w5 },
	{ w5, w4*w5, w5*w5 }
	};

	double Dx2[3][3] = { { 0,0,0 }  ,
	{ 0,1,w6 },
	{ 0, w6, w6*w6 } };

	double Dx3[3][3] = { { 0,0,0 }  ,
	{ 0,0,0 },
	{ 0,0,1 } };

	double Dx4[3][3] = { { 0, w4, 0 }  ,
	{ w1, 2 * w1*w4, w1*w5 },
	{ 0, w1*w5, 0 } };

	double Dx5[3][3] = { { 0,0,1 }  ,
	{ 0,0,w1*w4 },
	{ w1, w1*w4, 2 * w1*w5 } };

	double Dx6[3][3] = { { 0,0,0 }  ,
	{ 0,0, w2 },
	{ 0, w2, 2 * w2*w6 } };


	////Product of derivative of Diffusion tensor with Diffusion direction qk
	std::vector<std::vector<double>> qpartialDqsum;

	////Derivative of determinant of induced pullback metric gamma
	double detgamma_constant ;

	double sqrdetgamma = pow((detgamma_constant), -1 / 2);

	////Derivative of induced metric gamma
	std::vector<std::vector<std::vector<double>>> Derivgamma;

	///Derivative of determinant of induced metric
	std::vector<double> derivdetgamma;
	
	////Derivative of inverse induced metric gamma
	std::vector<std::vector<double> > DerivIIM;





												/////Member Functions

													//ELTensors

	//Embedded metric h_ij
	void hmatrixfunc(std::vector<double> &xhmatrix);   

	//Diffusion tensor product with diffusion direction q: qDq
	void qDqfunc(std::vector<double> &xqDqsum);

	///entries of pullback metric gamma
	void IMmatrixfunc(std::vector<std::vector<double> > &ximmatrix);


	////Induced Pullback metric gamma
	void detgammafunc(double &xdetgamma_constant);


	////Inverse of Induced metric gamma_nu,mu
	void IIMmatrixfunc(std::vector<std::vector<double> > &xinvimmatrix);

	///Christoffel Symbols
	void CSIComp(std::vector<std::vector<std::vector<double> > > &ChrisSym);



										//ELDerivTensors

		////Derivative of embedding map X
		void DerivXfunc(std::vector<std::vector<double>> &xDerivX);

		////Double derivative of embedding map X
		void DDerivXfunc(std::vector<std::vector<std::vector<double> > > &xDDerivX);


		////Derivative of induced metric
		void Derivgammafunc(std::vector<std::vector<std::vector<double>>> &xderivgamma);

		///Derivative of determinant of induced metric
		void Derivdetgammafunc(std::vector<double> &xderivdetgamma);

		////Derivative of inverse of induced metric; only for the terms needed in the EL equation
		void DerivIIMComp(std::vector<std::vector<double> > &derivIIM);


		/////Product of derivative of Diffusion tensor with Diffusion direction qk
		void qpartialDqfunc(std::vector<std::vector<double>> &xqpartialDqsum);


	
							////Euler Lagrange equations and Volume fraction

	////Initializing the tensors and their derivatives for the EL equations
	void TensorsandDerivTensorsInitialization();
	
	///Our EL eqns have three terms: 1)The sum involving the diffusion tensors , 
	const double term1Diff(int direction);

	/// 2)the partial derivatives of induced metric gamma with embedding map X and
	const double term2IX(int direction);

	///3) the term involving the Christoffel symbols.
	const double term3CS(int direction);

	///We will compute them separately and then add them for the EL scheme
	const double ELequation(int direction);

	//The iteration rule for the volume fraction
	const double VolfraIter();


};







#endif

