/** \file  ELGradient.cpp
\brief C++ source file for iterating the gradient descent scheme.

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
