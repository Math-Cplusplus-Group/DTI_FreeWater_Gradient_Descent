/** \file  CSVintomatrix.h
\brief C++ source file initializing matrices from csv.

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
*/#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

std::vector<std::vector<double> > CSVintoMatrix(std::ifstream& file, const int Graddirections, std::vector<std::vector<double> > &data, const int colsize = 3) {

	

	for (std::string line; std::getline(file, line); )
		//for (int row = 0; row != Graddirections; row++)
	{
		//std::string line;
		//std::getline(file, line);
		if (!file.good())
			break;

		std::stringstream iss(line);
		std::vector<double> rowset;
		std::string val;
		while (std::getline(iss, val, ','))
			//for (int col = 0; col != colsize; col++)
		{

			//std::string val;	
			//std::getline(iss, val, ',');
			//	if (!iss.good())
			//	break;

			//		std::stringstream convertor(val);
			//	convertor > > data[row][col];
			double s = std::stod(val);
			//std::cout << s;
			rowset.push_back(s);
			//	std::cout << s << ',';
		}

		//	std::cout << '\n' << rowset[0] << '\n';
		data.push_back(rowset);
		//std::cout << '\n' << data[0][1] << '\n';
	}



	return data;



}



