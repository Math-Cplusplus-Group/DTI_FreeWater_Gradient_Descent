/** \file  ELTensors.cpp
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

#include "ELTensors.h"
#include "ELInitialization.h"

////Embedded metric h_ij
void ELTensors::hmatrixfunc(std::vector<double> &xhmatrix) {
	
	//the induced metric h where the last entry is the cross term h_78=h_87
	xhmatrix.push_back(1);
	xhmatrix.push_back(1);
	xhmatrix.push_back(1);
       	xhmatrix.push_back(1 / (w1 * w1));
	xhmatrix.push_back(1 / (w2 * w2));
	xhmatrix.push_back(1 / (w3 * w3));
       	xhmatrix.push_back(2 * w1 * (w3 + w2 * w6 * w6) / (w2 * w3));
	xhmatrix.push_back(2 * w1 / w3);
      	xhmatrix.push_back(2 * w2 / w3);
	xhmatrix.push_back(-2 * w1 * w6 / w3);
	
}




//*********************************************************************************************************************************
//*********************************************************************************************************************************



//Product of Diffusin tensor with Diffusion direction qk
void ELTensors::qDqfunc(std::vector<double> &qDqsum) {
		
	qDqsum.resize(Eli.GradDirections);

	for (int k = 0; k != Eli.GradDirections; ++k) {
		for (int i = 0; i != 3; ++i)
		{
			for (int j = 0; j != 3; ++j)
			{
				qDqsum[k] += Eli(k,i) * Diff[i][j] * Eli(k, j);

			}
		}

	}

}



//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Induced metric gamma
void ELTensors::IMmatrixfunc(std::vector<std::vector<double> > &immatrix) {    

	immatrix.resize(3, std::vector<double>(3));

	for (int i = 0; i != 9; ++i)
	{
		//Wrt to DerivX[0]
		//the cross terms wrt to DerivX[0]
		immatrix[0][0] = 2 * DerivX[0][6] * DerivX[0][7] * hmatrix[8]; //the 2*(partial_1 X^7)*(partial_1 X^8)h^(78)
		immatrix[0][1] = (DerivX[0][6] * DerivX[1][7] + DerivX[0][7] * DerivX[1][6])* hmatrix[8];//the [(partial_1 X^8)*(partial_2 X^7)+(partial_1 X^7)*(partial_2 X^8)]*h^(78)
		immatrix[0][2] = (DerivX[0][6] * DerivX[2][7] + DerivX[0][7] * DerivX[2][6])* hmatrix[8];
		//std::cout << "This works5." << std::endl;

		//the diagonal terms
		immatrix[0][0] += DerivX[0][i] * DerivX[0][i] * hmatrix[i];
		immatrix[0][1] += DerivX[0][i] * DerivX[1][i] * hmatrix[i];
		immatrix[0][2] += DerivX[0][i] * DerivX[2][i] * hmatrix[i];
		//std::cout << "This works6." << std::endl;

		//wrt to DerivX[1]
		immatrix[1][0] = 2 * DerivX[1][6] * DerivX[1][7] * hmatrix[8]; //the 2*(partial_2 X^7)*(partial_2 X^8)h^(78)
		immatrix[1][2] = (DerivX[1][6] * DerivX[2][7] + DerivX[1][7] * DerivX[2][6])* hmatrix[8];


		//the diagonal terms
		immatrix[1][1] += DerivX[1][i] * DerivX[1][i] * hmatrix[i];
		immatrix[1][2] += DerivX[1][i] * DerivX[2][i] * hmatrix[i];

		//wrt to DerivX[2]
		immatrix[2][2] = 2 * DerivX[2][6] * DerivX[2][7] * hmatrix[8]; //the 2*(partial_2 X^7)*(partial_2 X^8)h^(78)

		immatrix[2][2] += DerivX[2][i] * DerivX[2][i] * hmatrix[i];
	}

	
}

//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Determinant of induced metric gamma

void ELTensors::detgammafunc(double &xdetgamma_constant_constant) {

	xdetgamma_constant_constant = IMmatrix[0][0] * (IMmatrix[1][1] * IMmatrix[2][2] - IMmatrix[1][2] * IMmatrix[1][2]) + IMmatrix[0][1] * (IMmatrix[0][2] * IMmatrix[1][2] - IMmatrix[0][1] * IMmatrix[2][2]) + IMmatrix[0][1] * (IMmatrix[0][1] * IMmatrix[1][2] - IMmatrix[0][2] * IMmatrix[1][1]);

}



//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Inverse of Induced metric


void ELTensors::IIMmatrixfunc(std::vector<std::vector<double> > &xInvIMmatrix) {
	
	xInvIMmatrix.resize(3, std::vector<double>(3));

	//entries of inverse induced metric gamma
	
	xInvIMmatrix[0][0] = (IMmatrix[1][1] * IMmatrix[2][2] - IMmatrix[1][2] * IMmatrix[1][2]) / detgamma_constant;
	xInvIMmatrix[1][1] = (IMmatrix[0][0] * IMmatrix[2][2] - IMmatrix[0][2] * IMmatrix[0][2]) / detgamma_constant;
	xInvIMmatrix[2][2] = (IMmatrix[0][0] * IMmatrix[1][1] - IMmatrix[0][1] * IMmatrix[0][1]) / detgamma_constant;

	xInvIMmatrix[0][1] = (IMmatrix[0][2] * IMmatrix[1][2] - IMmatrix[0][2] * IMmatrix[2][2]) / detgamma_constant;
	xInvIMmatrix[1][0] = xInvIMmatrix[0][1];

	xInvIMmatrix[0][2] = (IMmatrix[0][1] * IMmatrix[1][2] - IMmatrix[0][2] * IMmatrix[1][1]) / detgamma_constant;
	xInvIMmatrix[2][0] = xInvIMmatrix[0][2];


	xInvIMmatrix[1][2] = (IMmatrix[0][2] * IMmatrix[0][1] - IMmatrix[0][0] * IMmatrix[1][2]) / detgamma_constant;
	xInvIMmatrix[2][1] = xInvIMmatrix[1][2];

	
}



//*********************************************************************************************************************************
//*********************************************************************************************************************************

//Christoffel Symbols
//Gamma^i_{j,k}=1/2 h^(il)[partial_{j}h_{lk}+ partial_{k}h_{jl}-partial_{l}h_{jk}]

 void ELTensors::CSIComp(std::vector<std::vector<std::vector<double> > > &xChrisSym) {

	
	 xChrisSym.resize(9, std::vector<std::vector<double> >(9, std::vector<double>(9)));

	////For fixed i=4 the non-zero entries of the Gamma(4,j,k) matrix
	xChrisSym[3][3][3] = 1 / w1;
	xChrisSym[3][6][6] = -w1*w1*(w3 + w2*w6*w6) / (w2*w3);
	xChrisSym[3][6][7] = w1*w1*w6 / w3;
	xChrisSym[3][7][6] = w1*w1*w6 / w3;
	xChrisSym[3][7][7] = -w1*w1 / w3;

	////For fixed i=5
	xChrisSym[4][4][4] = 1 / w2;
	xChrisSym[4][6][6] = w1;
	xChrisSym[4][8][8] = -w2*w2 / w3;


	///For fixed i=6
	xChrisSym[5][5][5] = 1 / w3;
	xChrisSym[5][5][5] = w6*w6;
	xChrisSym[5][6][7] = -w1*w6;
	xChrisSym[5][7][6] = -w1*w6;
	xChrisSym[5][7][7] = w1;
	xChrisSym[5][8][8] = w2;


	///For fixed i=9
	xChrisSym[8][6][6] = -w1*w6 / w2;
	xChrisSym[8][6][7] = w1 / (2 * w2);
	xChrisSym[8][7][6] = w1 / (2 * w2);

}

