/** \file  ELTensors.cpp
\brief C++ source file initializing tensors.
Copyright 2016 by Andrew Colinet,Tomas Kojar

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
	ELTensors(std::vector<std::vector<std::vector<std::vector<double>>>> xCellX,
		double xdx,
		double xdy,
		double xdz,
		double xVolfn,
		ELInitialization xEli,
		std::vector<double> xAhat,
		std::ofstream& xmyfile
		) : CellX(xCellX), dx(xdx), dy(xdy), dz(xdz), Volfn(xVolfn), Eli(xEli), Ahat(xAhat),myfile(xmyfile) { };


	/// Copy constructor.



	//Member Functions

	//ELTensors

	//Embedded metric h_ij
	const std::vector<double> hmatrix();

	//Diffusion tensor product with diffusion direction q: qDq
	const double qDq(std::vector<double> qk);


	//Gamma^i_{j,k}=1/2 h^(il)[partial_{j}h_{lk}+ partial_{k}h_{jl}-partial_{l}h_{jk}]
	const double CSIComp(double i, double j, double k) const;

	//entries of gamma
	const std::vector<std::vector<double>> IMmatrix();

	//determinant gamma
	const double detgamma();

	//Induced metric gamma_nu,mu
	const std::vector<std::vector<double>> IIMComp();


	//Tensor derivatives

	//Derivative of embedding map X
	const std::vector<std::vector<double>> DerivX();

	//Second derivative of embedding map X
	const std::vector<std::vector<std::vector<double>>> DDerivX();

	//Derivative of induced metric gamma
	const std::vector<std::vector<std::vector<double>>> Derivgamma();

	//Deriv determinant gamma
	const double Derivdetgamma(double direction);


	//Derivative of inverse induced metric gamma_nu,mu
	const std::vector<std::vector<double>> DerivIIMComp();

	//Product of Derivative of Diffusion tensor D with kth diffusion direction q_k: (q_k^T) partialD (q_k)
	const double qpartialDq(double direction, std::vector<double> qk);


	//Euler Lagrange equations and Volume fraction

	///Our EL eqns have three terms: 1)The sum involving the diffusion tensors , 
	const double term1Diff(double direction);

	/// 2)the partial derivatives of induced metric gamma with embedding map X and
	const double term2IX(double direction);

	///3) the term involving the Christoffel symbols.
	const double term3CS(double direction);

	///We will compute them separately and then add them for the EL scheme
	const double ELequation(double direction);

	//The iteration rule for the volume fraction
	const double VolfraIter();


};







#endif

