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

#include "ELTensors.h"
#include "ELInitialization.h"

////Embedded metric h_ij
const std::vector<double> ELTensors::hmatrix() {
	//
	std::vector<double> X(9, 0);
	for (int i = 3; i != 9; i++) { X[i] = CellX[0][0][0][3]; }

	//the induced metric h where the last entry is the cross term h_78=h_87
	std::vector<double> hmetric = { 1,1,1, 1 / (X[3] * X[3]),1 / (X[4] * X[4]),1 / (X[5] * X[5]), 2 * X[3] * (X[5] + X[4] * X[8] * X[8]) / (X[4] * X[5]),2 * X[3] / X[5],  2 * X[4] / X[5],-2 * X[3] * X[8] / X[5] };

	return hmetric;

}



//*********************************************************************************************************************************
//*********************************************************************************************************************************



//Product of Diffusin tensor with Diffusion direction qk
const double ELTensors::qDq(std::vector<double> qk) {

	double w1 = CellX[0][0][0][3];
	double w2 = CellX[0][0][0][4];
	double w3 = CellX[0][0][0][5];
	double w4 = CellX[0][0][0][6];
	double 	w5 = CellX[0][0][0][7];
	double w6 = CellX[0][0][0][8];
	
	//std::cout << "This works3." << std::endl;

	double Diff[3][3] = { { w1, w1*w4, w1*w5 }  ,
	{ w1*w4, w2 + w1*w4*w4,w1*w4*w5 + w2*w6 },
	{ w1*w5, w1*w4*w5 + w2*w6, w3 + w1*w5*w5 + w2*w6*w6 } };

	double qDqsum = 0;

	for (int i = 0; i != 3; i++)
	{
		for (int j = 0; j != 3; j++)
		{
			qDqsum += qk[i] * Diff[i][j] * qk[j];
			//std::cout << "This works4." << std::endl;
		}
	}

	return qDqsum;
}



//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Induced metric gamma
const std::vector<std::vector<double>> ELTensors::IMmatrix() {

	std::vector<std::vector<double>> immatrix(3, std::vector<double>(3));

	//derivatives and image h_ij matrix	
	std::vector<double> partialxX(9, 0);
	std::vector<double> partialyX(9, 0);
	std::vector<double> partialzX(9, 0);
	std::vector<double> embeme(9, 0);

	for (int i = 0; i != 9; i++) {
		//std::cout << "Made it to for loop in IMmatrix" << " " << i << std::endl;
		partialxX[i] = (this->DerivX())[0][i];
		//std::cout << "This works." << std::endl;

		partialyX[i] = (this->DerivX())[1][i];
		//std::cout << "This works." << std::endl;

		partialzX[i] = (this->DerivX())[2][i];
		//std::cout << "This works." << std::endl;

		embeme[i] = (this->hmatrix())[i ];
		//std::cout << "This works now." << std::endl;

	}


	for (int i = 0; i != 9; i++)
	{
		//Wrt to partialxX
		//the cross terms wrt to partialxX
		immatrix[0][0] = 2 * partialxX[6] * partialxX[7] * embeme[8]; //the 2*(partial_1 X^7)*(partial_1 X^8)h^(78)
		immatrix[0][1] = (partialxX[6] * partialyX[7] + partialxX[7] * partialyX[6])* embeme[8];//the [(partial_1 X^8)*(partial_2 X^7)+(partial_1 X^7)*(partial_2 X^8)]*h^(78)
		immatrix[0][2] = (partialxX[6] * partialzX[7] + partialxX[7] * partialzX[6])* embeme[8];
		//std::cout << "This works5." << std::endl;

		//the diagonal terms
		immatrix[0][0] += partialxX[i] * partialxX[i] * embeme[i];
		immatrix[0][1] += partialxX[i] * partialyX[i] * embeme[i];
		immatrix[0][2] += partialxX[i] * partialzX[i] * embeme[i];
		//std::cout << "This works6." << std::endl;

		//wrt to partialyX
		immatrix[1][0] = 2 * partialyX[6] * partialyX[7] * embeme[8]; //the 2*(partial_2 X^7)*(partial_2 X^8)h^(78)
		immatrix[1][2] = (partialyX[6] * partialzX[7] + partialyX[7] * partialzX[6])* embeme[8];


		//the diagonal terms
		immatrix[1][1] += partialyX[i] * partialyX[i] * embeme[i];
		immatrix[1][2] += partialyX[i] * partialzX[i] * embeme[i];

		//wrt to partialzX
		immatrix[2][2] = 2 * partialzX[6] * partialzX[7] * embeme[8]; //the 2*(partial_2 X^7)*(partial_2 X^8)h^(78)

		immatrix[2][2] += partialzX[i] * partialzX[i] * embeme[i];
	}

	return immatrix;
}

//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Determinant of induced metric gamma

const double ELTensors::detgamma() {

	std::vector<std::vector<double>> IMComp(3,std::vector<double>(3));

	IMComp = this->IMmatrix();

	return IMComp[0][0] * (IMComp[1][1] * IMComp[2][2] - IMComp[1][2] * IMComp[1][2]) + IMComp[0][1] * (IMComp[0][2] * IMComp[1][2] - IMComp[0][1] * IMComp[2][2]) + IMComp[0][1] * (IMComp[0][1] * IMComp[1][2] - IMComp[0][2] * IMComp[1][1]);

}



//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Inverse of Induced metric


const std::vector<std::vector<double>> ELTensors::IIMComp() {
	double detgamma = (this->detgamma());
	std::vector<std::vector<double> > invimmatrix(3, std::vector<double>(3));

	std::vector<std::vector<double> > immatrix(3, std::vector<double>(3));

	immatrix = (this->IMmatrix());



	//entries of inverse induced metric gamma
	
	invimmatrix[0][0] = (immatrix[1][1] * immatrix[2][2] - immatrix[1][2] * immatrix[1][2]) / detgamma;
	invimmatrix[1][1] = (immatrix[0][0] * immatrix[2][2] - immatrix[0][2] * immatrix[0][2]) / detgamma;
	invimmatrix[2][2] = (immatrix[0][0] * immatrix[1][1] - immatrix[0][1] * immatrix[0][1]) / detgamma;

	invimmatrix[0][1] = (immatrix[0][2] * immatrix[1][2] - immatrix[0][2] * immatrix[2][2]) / detgamma;
	invimmatrix[1][0] = invimmatrix[0][1];

	invimmatrix[0][2] = (immatrix[0][1] * immatrix[1][2] - immatrix[0][2] * immatrix[1][1]) / detgamma;
	invimmatrix[2][0] = invimmatrix[0][2];


	invimmatrix[1][2] = (immatrix[0][2] * immatrix[0][1] - immatrix[0][0] * immatrix[1][2]) / detgamma;
	invimmatrix[2][1] = invimmatrix[1][2];

	return invimmatrix;
}



//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Christoffel Symbols
//Gamma^i_{j,k}=1/2 h^(il)[partial_{j}h_{lk}+ partial_{k}h_{jl}-partial_{l}h_{jk}]

const double ELTensors::CSIComp(double i, double j, double k) const {

	double w1 = CellX[0][0][0][3];
	double w2 = CellX[0][0][0][4];
	double w3 = CellX[0][0][0][5];
	double w6 = CellX[0][0][0][8];

	if (i == 4)
	{
		if (j == 4 && k == 4) return 1 / w1;

		else if (j == 7 && k == 7) return  -w1*w1*(w3 + w2*w6*w6) / (w2*w3);

		else if (j == 7 && k == 8) return  w1*w1*w6 / w3;

		else if (j == 8 && k == 7) return w1*w1*w6 / w3;
		else if (j == 8 && k == 8) return -w1*w1 / w3;

		else return 0;

	}
	else if (i == 5)
	{

		if (j == 5 && k == 5) return 1 / w2;

		else if (j == 7 && k == 7) return w1;

		else if (j == 9 && k == 9) return -w2*w2 / w3;

		else return 0;

	}
	else if (i == 6)
	{

		if (j == 6 && k == 6) return 1 / w3;

		else if (j == 7 && k == 7) return w6*w6;

		else if (j == 7 && k == 8) return	-w1*w6;
		else if (j == 8 && k == 7) return	-w1*w6;
		else if (j == 8 && k == 8) return w1;
		else if (j == 9 && k == 9) return	w2;


		else return 0;
	}
	else if (i == 9)
	{

		if (j == 7 && k == 7) return -w1*w6 / w2;

		else if (j == 7 && k == 8) return w1 / (2 * w2);

		else if (j == 8 && k == 7) return w1 / (2 * w2);

		else return 0;
	}

	else {

		return 0;
	}


}

