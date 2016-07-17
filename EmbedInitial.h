/** \file  EmbedInitial.h
\brief C++ header file initializing embedding map X.

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


#include "ELInitialization.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void Embedinitial(std::ofstream& myfile, std::vector<std::vector<double>> DiffGradDir, std::vector< std::vector<std::vector<std::vector<double>>> > Ahat, ELInitialization Eli, std::vector<std::vector<std::vector<double>>> volfn, std::vector< std::vector<std::vector<std::vector<double>>> > &embeddingmapXn){
		
	
	double Abitensor = 0;

	 embeddingmapXn.resize(9, std::vector<std::vector<std::vector<double> > >(Eli.nframesx, std::vector<std::vector<double>>(Eli.nframesy, std::vector<double >(Eli.nframesz))));
	
	 ////To use Least squares we express the 3x3 Diffusion matrix D for fixed x,y,z  as a 6-size vector (onlny 6 terms due to D's symmetry)
	
	////So DiffusionTensor={ DiffTenvec[0], DiffTenvec[1], DiffTenvec[2] ;  
	////					   DiffTenvec[1], DiffTenvec[3], DiffTenvec[4];
	////					   DiffTenvec[2], DiffTenvec[4], DiffTenvec[5] }



	//std::vector<std::vector<double>> Mcoeff(Eli.GradDirections, std::vector<double>(6));
	MatrixXf Mcoeff(Eli.GradDirections,6);
	VectorXf Ytiss(Eli.GradDirections);

	for (int x = 0; x != Eli.nframesx; x++) {
		for (int y = 0; y != Eli.nframesy; y++) {
			for (int z = 0; z != Eli.nframesz; z++) {

				for (int k = 0; k != Eli.GradDirections; ++k) {
					Mcoeff(k,0) = -Eli.bval*DiffGradDir[k][0] * DiffGradDir[k][0];
					Mcoeff(k,1) = -2 * Eli.bval * DiffGradDir[k][0] * DiffGradDir[k][1];
					Mcoeff(k,2) = -2 * Eli.bval * DiffGradDir[k][0] * DiffGradDir[k][2];
					Mcoeff(k,3) = -Eli.bval * DiffGradDir[k][1] * DiffGradDir[k][1];
					Mcoeff(k,4) = -2 * Eli.bval * DiffGradDir[k][1] * DiffGradDir[k][2];
					Mcoeff(k,5) = -2 * Eli.bval * DiffGradDir[k][2] * DiffGradDir[k][2];
					//myfile << '\n'<<"Mcoeff " << Mcoeff(k, 0) << ',' << Mcoeff(k, 1) << ',' << Mcoeff(k, 2) << ',' << Mcoeff(k, 3) << ',' << Mcoeff(k, 4) << ',' << Mcoeff(k, 5);
					Abitensor = (Ahat[x][y][z][k] - (1 - volfn[x][y][z])*Eli.Awater) /( volfn[x][y][z]+0.000001);

					Ytiss(k) = log(Abitensor+0.00000000001);
					//myfile << '\n' << ','<<"xyzk"<<','<<x<<y<<z<<k <<','<<"Ytiss" << ','<<Ytiss(k)<<',' << "volfn"<<','<< volfn[x][y][z] << ','<<"Ahat"<<','<< Ahat[x][y][z][k]<<',';
				}
				//std::cout << '\n' << "least squares Mcoeff" << Mcoeff;



											////Least squares to obtain the diffusion tensor: 
				
				//// Mcoeff*DiffTenvec= Ytiss => DiffTenvec=[ Mcoeff]^(-1)*Ytiss
				
				VectorXf DiffTenvec =Mcoeff.jacobiSvd(ComputeThinU | ComputeThinV).solve(Ytiss);
				
				////std::cout << '\n' << "least squares" << DiffTenvec << ',';
			
				////myfile <<"least squares" <<DiffTenvec(0) << ',' << DiffTenvec(1) << ',' << DiffTenvec(2)<< DiffTenvec(3)<< ',' << DiffTenvec(4)<< ',' << DiffTenvec(5)<<'\n';

										////Initializing the embedding map from the Diffusion tensor in Iwasawa coordinates

				////Identity map on the first three components
				embeddingmapXn[0][x][y][z] = x;
				embeddingmapXn[1][x][y][z] = y;
				embeddingmapXn[2][x][y][z] = z;

				//Fiber bundle components	
				//w1 term 
				embeddingmapXn[3][x][y][z] = DiffTenvec(0);  
				//w4 term
				embeddingmapXn[6][x][y][z] = DiffTenvec(1) / DiffTenvec(0);  
				//w5 term
				embeddingmapXn[7][x][y][z] = DiffTenvec(2) / DiffTenvec(0);  
				//w2 term
				embeddingmapXn[4][x][y][z] = DiffTenvec(3) - embeddingmapXn[3][x][y][z] * embeddingmapXn[6][x][y][z] * embeddingmapXn[6][x][y][z];  
				//w6 term
				embeddingmapXn[8][x][y][z] = (DiffTenvec(4) - embeddingmapXn[3][x][y][z] * embeddingmapXn[6][x][y][z] * embeddingmapXn[7][x][y][z]) / embeddingmapXn[4][x][y][z];
				//w3 term
				embeddingmapXn[5][x][y][z] = DiffTenvec(5) - embeddingmapXn[3][x][y][z] * embeddingmapXn[7][x][y][z] * embeddingmapXn[7][x][y][z] - embeddingmapXn[4][x][y][z] * embeddingmapXn[8][x][y][z] * embeddingmapXn[8][x][y][z];


	//std::cout << "embeddingmapXn " << embeddingmapXn[8][x][y][z] << '\n';
	//	myfile <<"embeddingmapXns: " <<embeddingmapXn[3][x][y][z] << ',' << embeddingmapXn[4][x][y][z] << ',' << embeddingmapXn[5][x][y][z] << ',' << embeddingmapXn[6][x][y][z] << ',' << embeddingmapXn[7][x][y][z] << ',' << embeddingmapXn[8][x][y][z] << '\n';
				
			}
		}
	}
	

}