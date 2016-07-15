/** \file  MainGradient.cpp
\brief C++ source file implementing classes for tensors.
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
#include "maxmintensor.h"
#include "CSVinto4Darray.h"
#include "CSVintoMatrix.h"
#include "CSVintoVector.h"
#include "MeanofTensor.h"
#include "Volumefraction.h"
#include "Volfnminmax.h"
#include "AhatInitializing.h"
#include "EmbedInitial.h"

int main() {

	std::ofstream myfile("LSresults.csv");
														//Physical variables
	int bval = 1000;
	double dwc = 0.003;
	//double Awater = exp(-bval* dwc);
	double Awater = 0;
	double lmax = 2.5 * 0.001; //constraint on maximum diffusion coefficient in tissue
	double lmin = 0.00001 ; //constraint on minimum diffusion coefficient in tissue
	double alpha = 1;


	//number of cells
	const int nuframesxtemp = 102;
	const int nuframesytemp = 102;
	const int nuframesx = 3;
	const int nuframesy = 3;
	const int nuframesz = 60;
	const int xGraddirections = 46;//Number of Gradient directions including the zero diffusion weight
	const int Graddirections = 44;//Without the zero diffusion weight
	//coordinates of voxel in deep White Matter (for initialization)
	//std::vector<int> WM { (int)(nuframesxtemp/2),(int)(nuframesytemp/2),(int)(nuframesz/2) };

	std::vector<int> WM{ 2,2,40};
	//coordinates of voxel with only free water(for initialization). Arbitrarily removed 3 points to be close to the boundary of brain.
	std::vector<int> WTR{ 2,2,50} ;

	std::cout << "line 39 "<<'\n';
	//Initializing DTI data
	std::ifstream mydtifile("DTIdataset3.csv");
	std::vector< std::vector<std::vector<std::vector<double>>> > Aatten = CSVinto4Darray(mydtifile, nuframesx, nuframesy, nuframesz, xGraddirections);
	std::cout << WM[0] << WM[1]<<WM[2]<< '\n';
	std::cout << "line 43 " << '\n';
	//Diffusion gradient Graddirections where the zero diffusion weight row vector was removed from the dataset
	std::ifstream mybvecfile("bvec2DTI30.csv");
	std::vector<std::vector<double>> DiffGradDir = CSVintoMatrix(mybvecfile, Graddirections);
	std::cout << "line 47" << '\n';
	//b-values vector in case they are not a single value but vary in each gradient direction
	//std::ifstream mybvalfile("bvalDTI30.csv");
	//std::vector<double> bval = CSVintoVector(mybvalfile, Graddirections);

	
	//initializing the physical data
	ELInitialization Eli(nuframesx, nuframesy, nuframesz, Graddirections, bval, DiffGradDir, alpha, Awater, lmin, lmax, WM, WTR);


	


													//Variables needed for iterations
											
													//mesh length
	double dt = 0.01;
	double dx = 0.01;
	double dy = 0.01;
	double dz = 0.01;
	double T = 5;
														
	//Attenuation vector at x,y,z with all Graddirections
	std::vector< std::vector<std::vector<std::vector<double>>> > Ahat=AhatInitializing(Aatten, nuframesx, nuframesy, nuframesz, Graddirections);

	std::cout << "line 72 " << '\n';
												//Volume Fraction initialization and min, max
	
	//mean over different diffusion Graddirections
	//std::vector<std::vector<std::vector<double>>> Amean = meanoftensor(Aatten, nuframesx, nuframesy, nuframesz, Graddirections);
		
	
	std::vector<std::vector<std::vector<std::vector<double>>>> volfnminmax = Volfnminmax(myfile,Ahat, Eli );
	std::vector<std::vector<std::vector<double>>> fmin = volfnminmax[0];
	std::vector<std::vector<std::vector<double>>> fmax = volfnminmax[1];
	std::cout << "line 82" << '\n';

	//Volume fraction at time t_n
	std::vector<std::vector<std::vector<double>>> Volfn = volumefraction(myfile,Aatten, Eli,fmin, fmax);
	
	//Volume fraction at time t_n+1
	std::vector<std::vector<std::vector<double>>> Volfnplus1 = Volfn;

																//Embedding map info

	//unit-Cell centered at X(x,y,z)
	std::vector<std::vector<std::vector<std::vector<double>>>> CellX(3, std::vector<std::vector<std::vector<double>>>(3, std::vector<std::vector<double>>(3, std::vector<double >(9, 0))));
	std::cout << "line 93" << '\n';
	
	//First coordinate is the X^i and the other three are x,y,z respectively for time t_n
	std::vector< std::vector<std::vector<std::vector<double>>> > embeddingmapXn;
	Embedinitial(myfile, DiffGradDir, Ahat, Eli, Volfn, embeddingmapXn);
	std::cout << "line 96" << '\n';
	//First coordinate is the X^i and the other three are x,y,z respectively for time t_n+1
	std::vector< std::vector<std::vector<std::vector<double>>> > embeddingmapXnplus1(9, std::vector<std::vector<std::vector<double>>>(nuframesx, std::vector<std::vector<double>>(nuframesy, std::vector<double >(nuframesz))));

	std::cout <<"least squares"<< embeddingmapXn[4][1][1][4];
	std::cout << "done";
	
	/*
	for (int i = 3; i != 9; i++) {
		for (int x = 0; x != nuframesx; x++) {
			for (int y = 0; y != nuframesy; y++) {
				for (int z = 0; z != nuframesz; z++) {

					myfile << embeddingmapXn[i][x][y][z] << ',';
				}
				myfile << '\n';
			}
			myfile << '\n';
		}
		myfile << '\n';
	}*/
	
		
																//Iteration
	for (int t = 0; t != T; t++) {
		for (int x = 1; x != nuframesx; x++) {
			for (int y = 1; y != nuframesy; y++) {
				for (int z = 1; z != nuframesz; z++) {

					//initialize the derivatives of X^i components. We need the x \pm dx and similarly for y,z.
					//std::cout << "this works "<< z;

					for (int i = 3; i != 9; i++) {

						for (int j1 = 0; j1 != 3; ++j1) {

							for (int j2 = 0; j2 != 3; ++j2) {

								for (int j3 = 0; j3 != 3; ++j3)
								{
									int m1 = j1;
									int m2 = j2;
									int m3 = j3;
									(j1 == 2) ? m1 = -1 : m1 = j1;
									(j2 == 2) ? m2 = -1 : m2 = j2;
									(j3 == 2) ? m3 = -1 : m3 = j3;
									if (z <= (nuframesz-2) && y<=(nuframesy - 2) && x<=(nuframesx-2)) {											
											CellX[j1][j2][j3][i] = embeddingmapXn[i][x + m1][y + m2][z + m3];	}

									else {
											if (x == (nuframesx - 1) && y != (nuframesy - 1) && z != (nuframesz - 1)) {
												CellX[j1][j2][j3][i] = embeddingmapXn[i][x][y + m2][z + m3];
											}


											else if (x != (nuframesx - 1) && y ==(nuframesy - 1) && z != (nuframesz - 1)) {
												CellX[j1][j2][j3][i] = embeddingmapXn[i][x + m1][y][z + m3];
											}

											else if (x != (nuframesx - 1) && y != (nuframesy - 1) && z == (nuframesz - 1)) {
												CellX[j1][j2][j3][i] = embeddingmapXn[i][x + m1][y+m1][z ];
											}

											else if (x == (nuframesx - 1) && y == (nuframesy - 1) && z != (nuframesz - 1)) {
												CellX[j1][j2][j3][i] = embeddingmapXn[i][x ][y][z+m3];

											}
											else if (x == (nuframesx - 1) && y != (nuframesy - 1) && z == (nuframesz - 1)) {
												CellX[j1][j2][j3][i] = embeddingmapXn[i][x ][y+m2][z];

											}

											else if (x != (nuframesx - 1) && y == (nuframesy - 1) && z == (nuframesz - 1)) {
												CellX[j1][j2][j3][i] = embeddingmapXn[i][x + m1][y][z];

											}

										}

									//std::cout << "this works for i,j1,j2,j3"<< i<< j1<< j2<< j3;

								}

							}

						}

					}



					std::cout << "this works for cell" ;

					for (int i = 3; i != 9; ++i) {

						

						ELTensors XandDX( CellX, dx, dy, dz, Volfn[x][y][z], Eli, Aatten[x][y][z] ,myfile);

					//	std::cout << "Object initialization worked." << std::endl;


						//run the finite difference calculation for embedding map X
						embeddingmapXnplus1[i][x][y][z] = embeddingmapXn[i][x][y][z] + dt*(XandDX.ELequation(i));
						//std::cout << "It worked." << std::endl;

						//run the finite difference calculation for volume fraction
						Volfnplus1[x][y][z] = Volfn[x][y][z] + dt*(XandDX.VolfraIter());
						//std::cout << Volfnplus1[x][y][z];
						//prepare for the next step
						embeddingmapXn[i][x][y][z] = embeddingmapXnplus1[i][x][y][z];
						Volfn[x][y][z] = std::min(std::max(Volfnplus1[x][y][z], fmin[x][y][z]), fmax[x][y][z]);

					}

					std::cout <<  ','<<"z is" <<','<< z<< ',' <<"y" <<y << ',' <<"x"<< x<<'\n';
				}

				
			}
			
		}

	
	}

	//Printing the Volume fraction and Diffustion tensor results into a CSV file
	std::ofstream myfile2("DTIresults.csv");
	//Creating diffusion tensor and printing it

	std::vector<std::vector< std::vector<std::vector<std::vector<double>>> >> DiffusionTensor;

	myfile << "(x,y,z)" << ',' << "Volume Fraction" << ',' << "Elements of Diffusion Tensor" << "\n";



	for (int x = 0; x != nuframesx; x++) {
		for (int y = 0; y != nuframesy; y++) {
			for (int z = 0; z != nuframesz; z++) {
				std::vector<std::vector<double>> DiffusionTensorxyz(3, std::vector<double>(3));

				DiffusionTensorxyz[0][0] = embeddingmapXnplus1[3][x][y][z];
				DiffusionTensorxyz[1][1] = embeddingmapXnplus1[4][x][y][z] + embeddingmapXnplus1[3][x][y][z] * embeddingmapXnplus1[6][x][y][z] * embeddingmapXnplus1[6][x][y][z];
				DiffusionTensorxyz[2][2] = embeddingmapXnplus1[5][x][y][z] + embeddingmapXnplus1[3][x][y][z] * embeddingmapXnplus1[7][x][y][z] * embeddingmapXnplus1[7][x][y][z] + embeddingmapXnplus1[4][x][y][z] * embeddingmapXnplus1[8][x][y][z] * embeddingmapXnplus1[8][x][y][z];

				DiffusionTensorxyz[0][1] = embeddingmapXnplus1[3][x][y][z] * embeddingmapXnplus1[6][x][y][z];
				DiffusionTensorxyz[2][1] = DiffusionTensorxyz[0][1];

				DiffusionTensorxyz[0][2] = embeddingmapXnplus1[3][x][y][z] * embeddingmapXnplus1[7][x][y][z];
				DiffusionTensorxyz[2][0] = DiffusionTensorxyz[0][2];

				DiffusionTensorxyz[1][2] = embeddingmapXnplus1[3][x][y][z] * embeddingmapXnplus1[6][x][y][z] * embeddingmapXnplus1[7][x][y][z] + embeddingmapXnplus1[4][x][y][z] * embeddingmapXnplus1[8][x][y][z];
				DiffusionTensorxyz[2][1] = DiffusionTensorxyz[1][2];


				myfile2 << "( " + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z) + ")" << Volfnplus1[x][y][z] << DiffusionTensorxyz[0][0] << ',' << DiffusionTensorxyz[0][1] << ',' << DiffusionTensorxyz[0][2] << ',' << DiffusionTensorxyz[1][0] << ',' << DiffusionTensorxyz[1][1] << ',' << DiffusionTensorxyz[1][2] << ',' << DiffusionTensorxyz[2][0] << ',' << DiffusionTensorxyz[2][1] << ',' << DiffusionTensorxyz[2][2] << "\n";

			}
		}
	}
	
	myfile.close();
	myfile2.close();
	
	
}
