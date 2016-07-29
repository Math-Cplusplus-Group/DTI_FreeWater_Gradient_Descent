/** \file  MeanofTensor.h
\brief C++ header file computing mean of 4Dtensors.

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

#include <vector>
#include <algorithm>

std::vector<std::vector<std::vector<double> > >  meanoftensor(std::vector< std::vector<std::vector<std::vector<double> > > > Aatten, const int nuframesx, const int nuframesy, const int nuframesz, double directions) {

	std::vector<std::vector<std::vector<double> > > result;
	double sum = 0;

	for (int x = 0; x != nuframesx; ++x) {
		std::vector<std::vector<double> > rowset2;
		for (int y = 0; y != nuframesy; ++y) {
			std::vector<double> rowset;
			for (int row = 0; row != nuframesz; ++row) {

				for (int k = 0; k != directions; ++k) {



					sum += Aatten[x][y][row][k] / directions;
					//std::cout << '\n'<<sum ;


				}
				rowset.push_back(sum);
			}
			rowset2.push_back(rowset);

		}
		result.push_back(rowset2);
	}


	return result;
}

