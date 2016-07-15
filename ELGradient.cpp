/** \file  ELGradient.cpp
\brief C++ source file for iterating the gradient descent scheme.
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


//Our EL eqns have three terms: 1)The sum involving the diffusion tensors , 
const double ELTensors::term1Diff(double compon) {
	double term1 = 0;
	double bval = Eli.bval;
	double  Awater = Eli.Awater;
	double alpha = Eli.alpha;
	double sqrdetgamma = pow((this->detgamma()), -1 / 2);

	//std::cout << "This works." << std::endl;

	for (int k = 0; k != Eli.Graddirections; k++) {
		std::vector<double> qk(3,0);
		for (int i = 0; i != qk.size(); i++) { qk[i] = Eli(k, i); };
		//std::cout << "This works1." << std::endl;
		double bqDq = bval*(this->qDq(qk));
		//std::cout << "qDq worked." << std::endl;

		term1 += (alpha*bval / sqrdetgamma) *(Volfn*exp(-bqDq) + (1 - Volfn)* Awater - (this->Ahat[k]))* exp(-bqDq)* (this->qpartialDq(compon, qk));

	}

	return term1;

}

// 2)the partial derivatives of induced metric gamma with embedding map X and
const double ELTensors::term2IX(double compon) {
	double term2 = 0;
	double sqrdetgamma = pow((this->detgamma()), -1 / 2);
	double detgamma = (this->detgamma());

	std::vector<double> partialX(3);
	for (int i = 0; i != 3; ++i) { partialX[i] = (this->DerivX())[i][compon];}
	
	std::vector<std::vector<double> >  partialinvgamma(3, std::vector<double>(3));
	partialinvgamma = (this->DerivIIMComp());
	
	std::cout << "43 ELG" << std::endl;
	std::vector<std::vector<double>> IIMCompmat(3, std::vector<double>(3));
	IIMCompmat = (this->IIMComp());

	std::cout << "47 ELG" << std::endl;
	std::vector<std::vector<double>> doubpartialX(3, std::vector<double>(3));
	doubpartialX = (this->DDerivX())[compon];

	std::cout << "50 ELG" << std::endl;
	for (int mu = 0; mu != 3; ++mu) {
		for (int nu = 0; nu != 3; ++nu) {
			std::cout << "53 ELG" << mu << nu << std::endl;
			term2 += -(1 / (2 * detgamma))*(this->Derivdetgamma(mu))*IIMCompmat[mu][nu]*partialX[nu] + partialinvgamma[mu][nu] *partialX[nu] + IIMCompmat[mu][nu] * doubpartialX[mu][nu];
			
		}
	}
	std::cout << "57 ELG" << std::endl;
	return term2;

}

//3) the term involving the Christoffel symbols.
const double ELTensors::term3CS(double direction) {

	double term3 = 0;

	std::vector<std::vector<double> > partialX(3, std::vector<double>(9));
	partialX = (this->DerivX());
	std::vector<std::vector<double> > IIMComp(3, std::vector<double>(3));
	IIMComp = (this->IIMComp());

	
	for (int mu = 0; mu != 3; mu++) {
		for (int nu = 0; nu != 3; nu++) {
			for (int j = 3; j != 9; j++) {
				for (int k = 3; k != 9; k++) {

					term3 += (this->CSIComp(direction, j, k))*IIMComp[mu][nu] * partialX[mu][j] * partialX[nu][direction];

					myfile << "term3"<< term3 <<','<<"IIMcomp"<< IIMComp[mu][nu] <<','<< partialX[mu][j]<< ','<< partialX[nu][direction] << '\n';
				}
			}
		}
	}

	return term3;
}

//We will compute them separately and then add them for the EL scheme
const double ELTensors::ELequation(double direction) {
	//std::cout << "Made it into ELequation" << std::endl;
	return term1Diff(direction) + term2IX(direction) + term3CS(direction);

}

//The iteration rule for the volume fraction
const double ELTensors::VolfraIter()
{
	double bval = Eli.bval;
	double Awater = Eli.Awater;
	double volfraiter = 0;
	for (int k = 0; k != Eli.Graddirections; k++)
	{
		std::vector<double> qk(3,0);
		for (int i = 0; i != qk.size(); i++) { qk[i] = Eli(k, i); };

		double bqDq = bval*(this->qDq(qk));

		volfraiter += -bval*(Volfn*exp(-bqDq) + (1 - Volfn)* Awater - Ahat[k])*(exp(-bqDq) - Awater);
	}

	return volfraiter;

}
