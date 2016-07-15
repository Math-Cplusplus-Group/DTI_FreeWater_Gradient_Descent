/** \file  CSVinto4Darray.h
\brief C++ header file initializing tensors from csv.
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
*/#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

std::vector< std::vector<std::vector<std::vector<double>>> > CSVinto4Darray(std::ifstream& file, const int nuframesx, const int nuframesy, const int nuframesz, const int Graddirections) {

	std::vector< std::vector<std::vector<std::vector<double> > > > data;

	for (int x = 0; x != nuframesx; ++x) {

		std::vector<std::vector<std::vector<double> > > rowset3;
		for (int y = 0; y != nuframesy; ++y) {

			std::vector<std::vector<double> > rowset2;
			for (int row = 0; row != nuframesz; ++row)
			{
				std::string line;
				std::getline(file, line);
				if (!file.good())
					break;

				std::stringstream iss(line);
				std::vector<double> rowset;
				for (int col = 0; col != Graddirections; ++col)
				{
					std::string val;
					std::getline(iss, val, ',');
					double s = std::stod(val);
					//if (!iss.good())
					//break;
					rowset.push_back(s);
					//			std::cout << s << ',';
					//std::stringstream convertor(val);
					//convertor >> data[x][y][row][col];

				}
				std::cout << '\n';
				rowset2.push_back(rowset);
				//std::cout << rowset2[0][2] << '\n';

			}
			rowset3.push_back(rowset2);
			//std::cout << rowset3[0][1][1] << '\n';
		}
		//	std::cout << '\n';
		data.push_back(rowset3);
		//std::cout << data[0][1][1][1]<<'\n';
	}

	return data;



}


