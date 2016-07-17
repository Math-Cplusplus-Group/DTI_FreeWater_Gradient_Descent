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

#include "ELTensors.h"
#include "ELInitialization.h"

//Derivative of embedding map X
void ELTensors::DerivXfunc(std::vector<std::vector<double>> &xDerivX) {

	xDerivX.resize(3, std::vector<double>(9));
	xDerivX[0][0] = 1; xDerivX[0][1] = 0; xDerivX[0][2] = 0;
	for (int i = 3; i != 9; ++i) {
		xDerivX[0][i] = (1 / dx)* (CellX[1][0][0][i] - CellX[0][0][0][i]);
	}

	xDerivX[1][0] = 0; xDerivX[1][1] = 1; xDerivX[1][2] = 0;
	for (int i = 3; i != 9; ++i) {
		xDerivX[1][i] = (1 / dx)* (CellX[0][1][0][i] - CellX[0][0][0][i]);
	}

	xDerivX[2][0] = 0; xDerivX[0][1] = 0; xDerivX[0][2] = 1;
	for (int i = 3; i != 9; ++i) {
		xDerivX[2][i] = (1 / dx)* (CellX[0][0][1][i] - CellX[0][0][0][i]);
	}
	

}

//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Double derivative of embedding map X

void ELTensors::DDerivXfunc(std::vector<std::vector<std::vector<double> > > &xDDerivX) {

	xDDerivX.resize(9, std::vector<std::vector<double> >(3, std::vector<double>  (3)));

	for (int i = 3; i != 9; ++i) {
		xDDerivX[i][0][0] = (1 / dx*dx)* (CellX[1][0][0][i] + CellX[2][0][0][i] - 2 * CellX[0][0][0][i]);
	
		xDDerivX[i][1][1] = (1 / dy*dy)* (CellX[0][1][0][i] + CellX[0][2][0][i] - 2 * CellX[0][0][0][i]);
	
		xDDerivX[i][2][2] = (1 / dz*dz)* (CellX[0][0][1][i] + CellX[0][0][2][i] - 2 * CellX[0][0][0][i]);
	
		xDDerivX[i][0][1] = (1 / dx*dy)* (CellX[1][1][0][i] - CellX[1][2][0][i] - CellX[2][1][0][i] + CellX[2][2][0][i]);
		xDDerivX[i][1][0] = xDDerivX[i][0][1];
		xDDerivX[i][0][2] = (1 / dx*dz)* (CellX[1][0][1][i] - CellX[1][0][2][i] - CellX[2][0][1][i] + CellX[2][0][2][i]);
		xDDerivX[i][2][0] = xDDerivX[i][0][2];
	
		xDDerivX[i][1][2] = (1 / dy*dz)* (CellX[0][1][1][i] - CellX[0][2][1][i] - CellX[0][1][2][i] + CellX[0][2][2][i]);
		xDDerivX[i][2][1] = xDDerivX[i][1][2];
	}
	

}


//********************************************************************************************************************************************
//********************************************************************************************************************************************

//Derivative of induced metric gamma
void ELTensors::Derivgammafunc(std::vector<std::vector<std::vector<double>>> &xDerivgamma) {
	
	xDerivgamma.resize(3, std::vector<std::vector<double>>(3, std::vector<double>(3)));

	std::vector<double> derivhmetric;

	for (int coord = 0; coord != 3; ++coord) {
		
		//derivative of induced metric
		derivhmetric = {
			0,
			0,
			0,
			-2 * DerivX[coord][3] / (w1 * w1 * w1),
			-2 * DerivX[coord][4] / (w2 * w2 * w2),
			-2 * DerivX[coord][5] / (w3 * w3 * w3),
			2 * (DerivX[coord][3] * (w3 + w2 * w6 * w6)*w2 * w3 + w1 * (DerivX[coord][5] + DerivX[coord][4] * w6 * w6 + 2 * w2 * w6 * DerivX[coord][8])*w2 * w3 - w1 * (w3 + w2 * w6 * w6)*(DerivX[coord][4] * w3 + w2 * DerivX[coord][5])) / (w2 * w2 * w3 * w3),
			2 * (DerivX[coord][3] * w3 - w1 * DerivX[coord][5]) / (w3 * w3),
			2 * (DerivX[coord][4] * w3 - w2 * DerivX[coord][5]) / (w3 * w3),
			-2 * (DerivX[coord][3] * w6 * w2 + w1 * DerivX[coord][8] * w2 - w2 * w6 * DerivX[coord][4]) / (w3 * w3)
		};



		for (int rowentries = 0; rowentries != 3; ++rowentries) {
			for (int colentries = 0; colentries != 3; ++colentries) {
				
				for (int i = 0; i != 9; ++i) {

					//the cross terms
					xDerivgamma[coord][rowentries][colentries] =	DDerivX[6][coord][rowentries] * DerivX[colentries][7] * hmatrix[8] + DerivX[rowentries][6] * DDerivX[7][coord][colentries] * hmatrix[8] + DerivX[rowentries][6] * DerivX[colentries][7] * derivhmetric[8] +

					DDerivX[7][coord][rowentries] * DerivX[colentries][6] * hmatrix[8] + DerivX[rowentries][7] * DDerivX[6][coord][colentries] * hmatrix[8] + DerivX[rowentries][7] * DerivX[colentries][6] * derivhmetric[8];

					//summing over the diagonal terms
					xDerivgamma[coord][rowentries][colentries] += DDerivX[i][coord][rowentries] * DerivX[colentries][i] * hmatrix[i] + DerivX[rowentries][i] * DDerivX[i][coord][colentries] * hmatrix[i] + DerivX[rowentries][i] * DerivX[colentries][i] * derivhmetric[i];
				}
			}
		}

	
	}



}


//Derivative of determinant of induced metric gamma
void ELTensors::Derivdetgammafunc(std::vector<double> &xderivdetgamma) {

	xderivdetgamma.resize(3);

	//computing the derivative of the deriminant using Jacobi's formula (detA)'=detA* tr( A^(-1) * (A)' )	
	for (int direction = 0; direction != 3; direction++) {
		for (int i1 = 0; i1 != 3; i1++) {
			for (int i2 = 0; i2 != 3; i2++)
			{
				xderivdetgamma[direction] += (detgamma_constant)* IMmatrix[i1][i2] * Derivgamma[direction][i2][i1];
			}
		}
	}
}




//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Derivative of inverse of induced metric gamma_nu,mu
void ELTensors::DerivIIMComp(std::vector<std::vector<double> > &xderivIIM) {


	xderivIIM.resize(3, std::vector<double>(3));

	//computing derivative of inverse using formula (A^-1)'=-(A^-1)* ( (A)' )* (A^-1)
		
	for (int coord = 0; coord != 3; ++coord) {
		for (int nu = 0; nu != 3; ++nu) {
			for (int i1 = 0; i1 != 3; ++i1) {
				for (int i2 = 0; i2 != 3; ++i2)
				{

					////std::cout << "line185 ELDer" << coord << " " << nu <<  i1<< i2<< std::endl;
					xderivIIM[coord][nu] -= InvIMmatrix[coord][i1] * Derivgamma[coord][i1][i2] * InvIMmatrix[i2][nu];

				}
			}
		}
	}




}


//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Product of derivative of Diffusion tensor with Diffusion direction qk
void ELTensors::qpartialDqfunc(std::vector<std::vector<double>> &xqpartialDqsum) {

	xqpartialDqsum.resize(Eli.GradDirections, std::vector<double>(9));

	for (int k = 0; k != Eli.GradDirections; k++) {
		for (int i = 0; i != 3; ++i)
		{
			for (int j = 0; j != 3; j++)
			{
				xqpartialDqsum[k][3] += Eli(k, i) * Dx1[i][j] * Eli(k, j);
				xqpartialDqsum[k][4] += Eli(k, i) * Dx2[i][j] * Eli(k, j);
				xqpartialDqsum[k][5] += Eli(k, i) * Dx3[i][j] * Eli(k, j);
				xqpartialDqsum[k][6] += Eli(k, i) * Dx4[i][j] * Eli(k, j);
				xqpartialDqsum[k][7] += Eli(k, i) * Dx5[i][j] * Eli(k, j);
				xqpartialDqsum[k][8] += Eli(k, i) * Dx6[i][j] * Eli(k, j);
			}
		}

	}
}