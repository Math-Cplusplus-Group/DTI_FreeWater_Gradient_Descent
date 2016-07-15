/** \file  Volfnminmax.h
C++ header file initializing the fmin and fmax of the free water volume fraction.
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
#include <vector>
#include <algorithm>
#include "ELInitialization.h"

std::vector< std::vector<std::vector<std::vector<double>>> >  Volfnminmax(std::ofstream& myfile, std::vector< std::vector<std::vector<std::vector<double>>> >  Ahat, ELInitialization Eli ) {
	double minAhat ;
	double maxAhat ;
	double bval = Eli.bval;
	double Awater = Eli.Awater;
	double lmin = Eli.lmin;
	double lmax = Eli.lmax;
	int nuframesx=Eli.nuframesx;
	int nuframesy= Eli.nuframesy;
	int nuframesz= Eli.nuframesz;
	int Graddirections=Eli.Graddirections;

	std::vector< std::vector<std::vector<std::vector<double>>> > result(2,std::vector<std::vector<std::vector<double>>>(nuframesx, std::vector<std::vector<double>>(nuframesy, std::vector<double >(nuframesz))));
	


	for (int x = 0; x != Eli.nuframesx; ++x) {

		for (int y = 0; y != Eli.nuframesy; ++y) {

			for (int z = 0; z != Eli.nuframesz; ++z) {
				
				//Computing minimum and maximum of Ahat for fixed x,y,z
				maxAhat = 0;
				minAhat = 1;
				for (int k = 0; k != Graddirections; ++k){
					double entry = Ahat[x][y][z][k];
					//std::cout << ',' << entry << ',' << k;
					//myfile << "entry" << ',' << entry << ','<< minAhat<<','<<maxAhat<<'\n';
					
					maxAhat = std::max(entry, maxAhat);
					minAhat = std::min(entry, minAhat);
					
					//myfile << "min and max "<<',' << minAhat << ',' << maxAhat << '\n';
					}

				result[0][x][y][z]= (minAhat - Awater) / (exp(-bval*lmin) - Awater);
				result[1][x][y][z] = (maxAhat - Awater) / (exp(-bval*lmax) - Awater);
				//myfile << "fmin" << ',' << result[0][x][y][z] << ','<< minAhat<<','<<Awater <<','<< exp(-bval*lmax)- Awater<< '\n';
				//myfile << "fmax" << ',' << result[1][x][y][z] << ',' << maxAhat - Awater << ',' << exp(-bval*lmin)<<','<< bval<< ',' << lmin<<','  << - Awater << '\n';

			}


		}

	}


	return result;
}

