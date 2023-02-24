//
//  testAIFMM.cpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#include "kernel.hpp"
#include "ACA.hpp"
#include "AIFMM.hpp"

int main(int argc, char* argv[]) {
	int sqrtRootN =	atoi(argv[1]);
	double L						=	atof(argv[2]); //domain size
	int TOL_POW         = atoi(argv[3]); // tolerance for compressions in powers of 10
	int nParticlesInLeafAlong1D	=	atoi(argv[4]); // assuming the particles are located at tensor product chebyNodes/uniform

	int nLevelsUniform		=	ceil(log(double(sqrtRootN)/nParticlesInLeafAlong1D)/log(2));

	std::vector<pts2D> particles_X, particles_Y; //dummy variables
	userkernel* mykernel		=	new userkernel(sqrtRootN, L, nLevelsUniform);
	unsigned N = sqrtRootN*sqrtRootN;
	unsigned nLevels            = mykernel->nLevelsResolve;//for quad tree

	double start, end;
	int n_Dimension = 2;
	double* locations       =       new double [N*n_Dimension];   //      Stores all the locations.

	std::cout << std::endl << "N: " << N << std::endl;

	std::cout << std::endl << "nLevelsUniform: " << nLevelsUniform << std::endl;

	// Generates random locations and random values of property. Change the locations and rhs as per your need
	/////////////////////////////// GIVE INPUTS HERE ///////////////////////////////////////////////
	unsigned count_Location =       0;
	for (unsigned j=0; j<N; ++j) {
			for (unsigned k=0; k<n_Dimension; ++k) {
				if(k == 0) {
					locations[count_Location]       =       mykernel->gridPoints[j].x;//2*(int(rand())%2)-1;//double(rand())/double(RAND_MAX);
					++count_Location;
				}
				else {
					locations[count_Location]       =       mykernel->gridPoints[j].y;//2*(int(rand())%2)-1;//double(rand())/double(RAND_MAX);
					++count_Location;
				}
			}
	}
	Vec aifmm_x;
	double err;
	//////////// ONES AT RANDOM PLACES  //////////////////////////
	// For testing the code, x is considered to be vector of ones and zeros and the corresponding b is calculated in the following code snippet
	{
		Vec true_x=Vec::Zero(N);
		int n = N/500; //randomly choosing n different indices where b is set to 1, b at the rest of the indices is set to 0
		srand(time(NULL));
		std::set<int> s;
		while(s.size() < n) {
			int index	=	rand()%N;
			s.insert(index);
		}
		std::set<int>::iterator it;
		for (it = s.begin(); it != s.end(); it++) {
			true_x(*it) = 1.0;
		}
		Vec true_Ax = Vec::Zero(N);
		for (it = s.begin(); it != s.end(); it++) {
			true_Ax = true_Ax + mykernel->getCol(N,*it);
		}
	}
	int TOL_POW_NCA_AIFMM = TOL_POW;
	start	=	omp_get_wtime();
	AIFMM<userkernel> *aifmm = new AIFMM<userkernel>(mykernel, N, nLevelsUniform, TOL_POW_NCA_AIFMM, locations);
	end		=	omp_get_wtime();
	double timeAssemble =	(end-start);

	//////////////////////// ELIMINATE PHASE ////////////////////////////////////////
	start	=	omp_get_wtime();
	aifmm->factorize();
	end		=	omp_get_wtime();
	double timeFactorize =	(end-start);

	aifmm->backSubstitute1(true_Ax); //assigning rhs

	start	=	omp_get_wtime();
	aifmm->backSubstitute2();// doing schur compliment operations on rhs (related to factorization)
	end		=	omp_get_wtime();
	timeFactorize +=	(end-start);

	start	=	omp_get_wtime();
	aifmm->backSubstitute3();// back substitution phase or the solve phase
	end		=	omp_get_wtime();
	double timeSolve =	(end-start);

	aifmm->getSolution(aifmm_x);

	err = (true_x - aifmm_x).norm()/true_x.norm();

	std::cout << std::endl << "Max rank of compressible sub-matrices: " << aifmm->A->getMaxRank() << std::endl;  // (including those of compressible fill-ins)

	std::cout << std::endl << "Assemble time: " << timeAssemble << std::endl;

	std::cout << std::endl << "Factorization time: " << timeFactorize << std::endl;

	std::cout << std::endl << "Solve time: " << timeSolve << std::endl;

	std::cout << std::endl << "Relative forward error in the solution: " << err << std::endl; // relative forward error in 2 norm sense

	delete aifmm;

	delete locations;
	delete mykernel;
}
