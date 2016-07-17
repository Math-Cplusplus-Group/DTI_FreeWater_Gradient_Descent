/** \file   Volumefraction.h
C++ header file initializing the free water volume fraction.
Copyright 2016 by Tomas Kojar

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
std::vector<std::vector<std::vector<double>>>  volumefraction(std::ofstream& myfile, std::vector< std::vector<std::vector<std::vector<double>>> > Aatten, ELInitialization Eli, std::vector<std::vector<std::vector<double>>> fmin, std::vector<std::vector<std::vector<double>>> fmax) {
	
	int nuframesx = Eli.nframesx;
	int nuframesy = Eli.nframesy;
	int nuframesz = Eli.nframesz;
	int Graddirections = Eli.GradDirections;
	std::vector<int> WM = Eli.WM;
	std::vector<int> WTR = Eli.WTR;
	
	//Stis is a baseline value from a voxel expected not to have a free water compartment(typically within deep white matter structures)
	//Swat is an intensity value from a voxel containing only a free water component
	double Stis=Aatten[WM[0]][WM[1]][WM[2]][0];
	double Swat = Aatten[WTR[0]][WTR[1]][WTR[2]][0];
	//myfile << Swat << ',' << Stis << '\n';

	std::vector<std::vector<std::vector<double>>> result(std::vector<std::vector<std::vector<double>>>(nuframesx, std::vector<std::vector<double>>(nuframesy, std::vector<double >(nuframesz))));
	
	double s=0;

	for (int x = 0; x != nuframesx; ++x) {
		
		for (int y = 0; y != nuframesy; ++y) {

			for (int z = 0; z != nuframesz; ++z) {

				s= 1 - log(Aatten[x][y][z][0] / Stis) / log(Swat / Stis);
				myfile << '\n' << "Volfn before" << ',' << s;
				
				if (s<=fmin[x][y][z] || s>=fmax[x][y][z]) { result[x][y][z] = (fmin[x][y][z] + fmin[x][y][z])*0.5; }
				else{ result[x][y][z] = 1 - log(Aatten[x][y][z][0] / Stis) / log(Swat / Stis); }

				myfile << '\n' << "Volfn" << ','<<result[x][y][z] <<"fmin"<<','<< fmin[x][y][z]<<','<< "fmax" << ',' << fmax[x][y][z]<<','<<"xyz"<<','<<x<<y<<z;
								
			}
		}
	}


	return result;
}

