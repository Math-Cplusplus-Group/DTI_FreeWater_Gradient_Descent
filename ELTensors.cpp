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

//*********************************************************************************************************************************
//*********************************************************************************************************************************

//*********************************************************************************************************************************
//*********************************************************************************************************************************

//*********************************************************************************************************************************
//*********************************************************************************************************************************
////Initializing the tensors and their derivatives for the EL equations

void ELTensors::TensorsandDerivTensorsInitialization(){

	///The order of calling the functions matters, since some of them depend on others


								///Foundational Tensors not depending on other tensors
	//Embedded metric h_ij
	hmatrixfunc( hmatrix);

	///Christoffel Symbols
	CSIComp(ChrisSym);

	//Diffusion tensor product with diffusion direction q: qDq
	qDqfunc(qDqsum);
	
	/////Product of derivative of Diffusion tensor with Diffusion direction qk
	qpartialDqfunc(qpartialDqsum);

	////Derivative of embedding map X
	DerivXfunc(DerivX);

	////Double derivative of embedding map X
	DDerivXfunc(DDerivX);

	

											////Tensors that depend on other tensors
	
	///entries of pullback metric gamma. It depends on DerivX.
	IMmatrixfunc( IMmatrix);


	////Induced Pullback metric gamma. It depends on IMmatrix.
	detgammafunc( detgamma_constant);

	
	////Inverse of Induced metric gamma_nu,mu. It depends on IMmatrix and detgamma_constant.
	 IIMmatrixfunc( InvIMmatrix);

	 
	 ////Derivative of induced metric. It depends on DerivX, DDerivX, hmatrix.
	 Derivgammafunc(Derivgamma);
	
	 ///Derivative of determinant of induced metric. It depends on IMmatrix, detgamma_constant, Derivgamma.
	 Derivdetgammafunc(derivdetgamma);

	 ////Derivative of inverse of induced metric; only for the terms needed in the EL equation. It depends on InvIMmatrix, Derivgamma.
	 DerivIIMComp(DerivIIM);






}






//Our EL eqns have three terms: 1)The sum involving the diffusion tensors , 
const double ELTensors::term1Diff(int compon) {
	double term1 = 0;
	double bqDq = 0 ;

	for (int k = 0; k != Eli.GradDirections; k++) {
		bqDq = bval*qDqsum[k];
		term1 += (alpha*bval / sqrdetgamma) *(Volfn*exp(-bqDq) + (1 - Volfn)* Awater - (this->Ahat[k]))* exp(-bqDq)* qpartialDqsum[k][compon];
	}

	return term1;

}

// 2)the partial derivatives of induced metric gamma with embedding map X and
const double ELTensors::term2IX(int compon) {
	double term2 = 0;
	
	
	for (int mu = 0; mu != 3; ++mu) {
		for (int nu = 0; nu != 3; ++nu) {
		////std::cout << "53 ELG" << mu << nu << std::endl;
			term2 += -(1 / (2 * detgamma_constant))*derivdetgamma[mu] * InvIMmatrix[mu][nu]*DerivX[nu][compon] + DerivIIM[mu][nu] * DerivX[nu][compon] + InvIMmatrix[mu][nu] * DDerivX[compon][mu][nu];
			
		}
	}

	return term2;

}

//3) the term involving the Christoffel symbols.
const double ELTensors::term3CS(int direction) {

	double term3 = 0;
		
	for (int mu = 0; mu != 3; mu++) {
		for (int nu = 0; nu != 3; nu++) {
			for (int j = 3; j != 9; j++) {
				for (int k = 3; k != 9; k++) {

					term3 += ChrisSym[direction][j][k]*InvIMmatrix[mu][nu] * DerivX[mu][j] * DerivX[nu][direction];

	///myfile << "term3"<< term3 <<','<<"IIMcomp"<< IIMComp[mu][nu] <<','<< partialX[mu][j]<< ','<< partialX[nu][direction] << '\n';
				}
			}
		}
	}

	return term3;
}


//We will compute them separately and then add them for the EL scheme
const double ELTensors::ELequation(int direction) {
	
	return term1Diff(direction)+ term2IX(direction) + term3CS(direction);

}

//The iteration rule for the volume fraction
const double ELTensors::VolfraIter()
{	
	double volfraiter = 0;
	double bqDq = 0;
	
	for (int k = 0; k != Eli.GradDirections; k++)
	{	
		 bqDq = bval*qDqsum[k];
		volfraiter += -bval*(Volfn*exp(-bqDq) + (1 - Volfn)* Awater - Ahat[k])*(exp(-bqDq) - Awater);
	}

	return volfraiter;

}



//*********************************************************************************************************************************
//*********************************************************************************************************************************

//*********************************************************************************************************************************
//*********************************************************************************************************************************

//*********************************************************************************************************************************
//*********************************************************************************************************************************



//Derivative of embedding map X
void ELTensors::DerivXfunc(std::vector<std::vector<double> > &xDerivX) {

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
void ELTensors::Derivgammafunc(std::vector<std::vector<std::vector<double> > > &xDerivgamma) {
	
	xDerivgamma.resize(3, std::vector<std::vector<double> >(3, std::vector<double>(3)));

	std::vector<double> derivhmetric;

	for (int coord = 0; coord != 3; ++coord) {
		
		//derivative of induced metric
		derivhmetric.push_back(0);
		derivhmetric.push_back(0);
		derivhmetric.push_back(0);
		derivhmetric.push_back(-2 * DerivX[coord][3] / (w1 * w1 * w1));
		derivhmetric.push_back(-2 * DerivX[coord][4] / (w2 * w2 * w2));
		derivhmetric.push_back(-2 * DerivX[coord][5] / (w3 * w3 * w3));
		derivhmetric.push_back(2 * (DerivX[coord][3] * (w3 + w2 * w6 * w6)*w2 * w3 + w1 * (DerivX[coord][5] + DerivX[coord][4] * w6 * w6 + 2 * w2 * w6 * DerivX[coord][8])*w2 * w3 - w1 * (w3 + w2 * w6 * w6)*(DerivX[coord][4] * w3 + w2 * DerivX[coord][5])) / (w2 * w2 * w3 * w3));
		derivhmetric.push_back(2 * (DerivX[coord][3] * w3 - w1 * DerivX[coord][5]) / (w3 * w3));
		derivhmetric.push_back(2 * (DerivX[coord][4] * w3 - w2 * DerivX[coord][5]) / (w3 * w3));
		derivhmetric.push_back(-2 * (DerivX[coord][3] * w6 * w2 + w1 * DerivX[coord][8] * w2 - w2 * w6 * DerivX[coord][4]) / (w3 * w3));



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
void ELTensors::qpartialDqfunc(std::vector<std::vector<double> > &xqpartialDqsum) {

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



