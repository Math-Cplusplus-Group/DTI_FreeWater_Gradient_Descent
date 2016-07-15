/** \file  ELDerivTensors.cpp
\brief C++ source file initializing the derivatives of tensors.
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

#include <cmath>
#include <vector>
#include "ELTensors.h"
#include "ELInitialization.h"


//Derivative of embedding map X
const std::vector<std::vector<double> > ELTensors::DerivX() {

	std::vector<std::vector<double>> derivX(3, std::vector<double>(9));

	derivX[0][0] = 1; derivX[0][1] = 0; derivX[0][2] = 0;
	for (int i = 3; i != 9; ++i) {
		derivX[0][i] = (1 / dx)* (CellX[1][0][0][i] - CellX[0][0][0][i]);
	}

	derivX[1][0] = 0; derivX[1][1] = 1; derivX[1][2] = 0;
	for (int i = 3; i != 9; ++i) {
		derivX[1][i] = (1 / dx)* (CellX[0][1][0][i] - CellX[0][0][0][i]);
	}

	derivX[2][0] = 0; derivX[0][1] = 0; derivX[0][2] = 1;
	for (int i = 3; i != 9; ++i) {
		derivX[2][i] = (1 / dx)* (CellX[0][0][1][i] - CellX[0][0][0][i]);
	}


	return derivX;

}

//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Double derivative of embedding map X

const std::vector<std::vector<std::vector<double>>> ELTensors::DDerivX() {

	std::vector<std::vector<std::vector<double> > > dderiv(9, std::vector<std::vector<double>>(3, std::vector<double>(3)));
	
	for (int i = 3; i != 9; ++i) {
		dderiv[i][0][0] = (1 / dx*dx)* (CellX[1][0][0][i] + CellX[2][0][0][i] - 2 * CellX[0][0][0][i]);
	}

	

	for (int i = 3; i != 9; ++i) {
		dderiv[i][1][1] = (1 / dy*dy)* (CellX[0][1][0][i] + CellX[0][2][0][i] - 2 * CellX[0][0][0][i]);
	}

	



	for (int i = 3; i != 9; ++i) {
		dderiv[i][2][2] = (1 / dz*dz)* (CellX[0][0][1][i] + CellX[0][0][2][i] - 2 * CellX[0][0][0][i]);
	}


	for (int i = 3; i != 9; ++i) {
		dderiv[i][0][1] = (1 / dx*dy)* (CellX[1][1][0][i] - CellX[1][2][0][i] - CellX[2][1][0][i] + CellX[2][2][0][i]);
		dderiv[i][1][0] = dderiv[i][0][1];
	}

	for (int i = 3; i != 9; ++i) {
		dderiv[i][0][2] = (1 / dx*dz)* (CellX[1][0][1][i] - CellX[1][0][2][i] - CellX[2][0][1][i] + CellX[2][0][2][i]);
		dderiv[i][2][0] = dderiv[i][0][2];
	}
	
	for (int i = 3; i !=9; ++i) {
		dderiv[i][1][2] = (1 / dy*dz)* (CellX[0][1][1][i] - CellX[0][2][1][i] - CellX[0][1][2][i] + CellX[0][2][2][i]);
		dderiv[i][2][1] = dderiv[i][1][2];
	}



	return dderiv;

}


//********************************************************************************************************************************************
//********************************************************************************************************************************************

//Derivative of induced metric gamma
const std::vector<std::vector<std::vector<double>>> ELTensors::Derivgamma() {

	std::vector<std::vector<std::vector<double> > > derivgamma(9, std::vector<std::vector<double>>(3, std::vector<double>(3)));

	//the embedding map where we start with X[0]=0 and X[i]=X^i for readability
	std::vector<double> X(9, 0);
	for (int i = 0; i !=9; ++i) { X[i] = CellX[0][0][0][i]; }

	//image metric h
	std::vector<double> h = this->hmatrix();

	//derivative of X: eg. co=0 means partial X /partial x
	std::vector<std::vector<double> > derivX(3, std::vector<double>(9));
	for (int co = 0; co != 3; ++co) {
		for (int i = 0; i != 9; ++i) { derivX[co][i] = (this->DerivX())[co][i]; }
	}

	
	for (int coord = 0; coord != 3; ++coord) {
		//derivative of induced metric
		std::vector<double> derivhmetric = {
			0,
			0,
			0,
			-2 * derivX[coord][3] / (X[3] * X[3] * X[3]),
			-2 * derivX[coord][4] / (X[4] * X[4] * X[4]),
			-2 * derivX[coord][5] / (X[5] * X[5] * X[5]),
			2 * (derivX[coord][3] * (X[5] + X[4] * X[8] * X[8])*X[4] * X[5] + X[3] * (derivX[coord][5] + derivX[coord][4] * X[8] * X[8] + 2 * X[4] * X[8] * derivX[coord][8])*X[4] * X[5] - X[3] * (X[5] + X[4] * X[8] * X[8])*(derivX[coord][4] * X[5] + X[4] * derivX[coord][5])) / (X[4] * X[4] * X[5] * X[5]),
			2 * (derivX[coord][3] * X[5] - X[3] * derivX[coord][5]) / (X[5] * X[5]),
			2 * (derivX[coord][4] * X[5] - X[4] * derivX[coord][5]) / (X[5] * X[5]),
			-2 * (derivX[coord][3] * X[8] * X[4] + X[3] * derivX[coord][8] * X[4] - X[4] * X[8] * derivX[coord][4]) / (X[5] * X[5])
		};
		


		for (int rowentries = 0; rowentries != 3; ++rowentries) {
			for (int colentries = 0; colentries != 3; ++colentries) {
				std::vector<std::vector<std::vector<double>>> dderivX(9, std::vector<std::vector<double>>(3, std::vector<double>(3)));
				dderivX = this->DDerivX();

				for (int i = 0; i != 9; ++i) {

					//the cross terms
					derivgamma[coord][rowentries][colentries] =

						dderivX[6][coord][rowentries] * derivX[colentries][7] * h[8] + derivX[rowentries][6] * dderivX[7][coord][colentries] * h[8] + derivX[rowentries][6] * derivX[colentries][7] * derivhmetric[8] +

						dderivX[7][coord][rowentries] * derivX[colentries][6] * h[8] + derivX[rowentries][7] * dderivX[6][coord][colentries] * h[8] + derivX[rowentries][7] * derivX[colentries][6] * derivhmetric[8];

					//summing over the diagonal terms
					derivgamma[coord][rowentries][colentries] += dderivX[i][coord][rowentries] * derivX[colentries][i ] * h[i] + derivX[rowentries][i ] * dderivX[i][coord][colentries] * h[i] + derivX[rowentries][i] * derivX[colentries][i ] * derivhmetric[i];
				}
			}
		}

		//std::cout << "41 ELD coord= " << coord << std::endl;


	}

	return derivgamma;

}


//Derivative of determinant of induced metric gamma
const double ELTensors::Derivdetgamma(double direction) {
	double derivdet = 0;

	//computing the derivative of the deriminant using Jacobi's formula (detA)'=detA* tr( A^(-1) * (A)' )

	for (int i1 = 0; i1 != 3; i1++) {
		for (int i2 = 0; i2 != 3; i2++)
		{
			derivdet += (this->detgamma()) * (this->IIMComp())[i1][i2] * (this->Derivgamma())[direction][i2][i1];
		}
	}
	return derivdet;
}




//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Derivative of inverse of induced metric gamma_nu,mu
const std::vector<std::vector<double>> ELTensors::DerivIIMComp() {

	std::vector<std::vector<double>> derivIIM(3, std::vector<double>(3));
	//computing derivative of inverse using formula (A^-1)'=-(A^-1)* ( (A)' )* (A^-1)
//	std::cout << "New Message1" << std::endl;

	std::vector<std::vector<std::vector<double> > > derivgamma = (this->Derivgamma());
	std::cout << "line180 ";
	for (int coord = 0; coord != 3; ++coord) {
		for (int nu = 0; nu != 3; ++nu) {
			for (int i1 = 0; i1 != 3; ++i1) {
				for (int i2 = 0; i2 != 3; ++i2)
				{

					std::cout << "line185 " << coord << " " << nu <<  std::endl;
					derivIIM[coord][nu] -= (this->IIMComp())[coord][i1]*  derivgamma[coord][i1][i2]*(this->IIMComp())[i2][nu];

				}
			}
		}
	}


	return derivIIM;

}


//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Product of derivative of Diffusion tensor with Diffusion direction qk
const double ELTensors::qpartialDq(double direction, std::vector<double> qk) {

	double w1 = CellX[0][0][0][3];
	double w2 = CellX[0][0][0][4];
	double w3 = CellX[0][0][0][5];
	double w4 = CellX[0][0][0][6];
	double 	w5 = CellX[0][0][0][7];
	double w6 = CellX[0][0][0][8];

	if (direction == 4) {

		double Dx1[3][3] = { { 1, w4, w5 }  ,
		{ w4, w4*w4, w4*w5 },
		{ w5, w4*w5, w5*w5 }
		};

		double qDqsum = 0;

		for (int i = 0; i != 3; ++i)
		{
			for (int j = 0; j != 3; j++)
			{
				qDqsum += qk[i] * Dx1[i][j] * qk[j];
			}
		}

		return qDqsum;
	}


	else if (direction == 5) {

		double Dx2[3][3] = { { 0,0,0 }  ,
		{ 0,1,w6 },
		{ 0, w6, w6*w6 } };

		double qDqsum = 0;

		for (int i = 0; i != 3; ++i)
		{
			for (int j = 0; j != 3; j++)
			{
				qDqsum += qk[i] * Dx2[i][j] * qk[j];
			}
		}

		return qDqsum;

	}


	else if (direction == 6) {

		double Dx3[3][3] = { { 0,0,0 }  ,
		{ 0,0,0 },
		{ 0,0,1 } };

		double qDqsum = 0;

		for (int i = 0; i != 3; ++i)
		{
			for (int j = 0; j != 3; j++)
			{
				qDqsum += qk[i] * Dx3[i][j] * qk[j];
			}
		}
		return qDqsum;
	}


	else if (direction == 7) {

		double Dx4[3][3] = { { 0, w4, 0 }  ,
		{ w1, 2 * w1*w4, w1*w5 },
		{ 0, w1*w5, 0 } };
		double qDqsum = 0;

		for (int i = 0; i != 3; ++i)
		{
			for (int j = 0; j != 3; j++)
			{
				qDqsum += qk[i] * Dx4[i][j] * qk[j];
			}
		}

		return qDqsum;
	}

	else if (direction == 8) {

		double Dx5[3][3] = { { 0,0,1 }  ,
		{ 0,0,w1*w4 },
		{ w1, w1*w4, 2 * w1*w5 } };
		double qDqsum = 0;

		for (int i = 0; i != 3; ++i)
		{
			for (int j = 0; j != 3; j++)
			{
				qDqsum += qk[i] * Dx5[i][j] * qk[j];
			}
		}
		return qDqsum;

	}

	else if (direction == 9) {

		double Dx6[3][3] = { { 0,0,0 }  ,
		{ 0,0, w2 },
		{ 0, w2, 2 * w2*w6 } };
		double qDqsum = 0;

		for (int i = 0; i != 3; ++i)
		{
			for (int j = 0; j != 3; j++)
			{
				qDqsum += qk[i] * Dx6[i][j] * qk[j];
			}
		}
		return qDqsum;

	}
	else return 0;
}
