//
//  testFMM2D.cpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#include "kernel.hpp"
#include "ACA.hpp"
#include "AIFMM.hpp"
#include "AFMM.hpp"
#include <filesystem>
#include "gmres.hpp"

int main(int argc, char* argv[]) {
	#ifdef yestoLS
		int nChebNodes			=	atoi(argv[1]);//number of gridPoints in leaf box in 1D
		int treeAdaptivity	=	1;//tolerance for tree adpativity
		kappa 							= atof(argv[2]);//frequency of operation
		int degreeOfBases 	= atoi(argv[3]); //for discretization of Lippmann-Schwinger equation
		double L						=	atof(argv[4]); //domain size
		Qchoice		    			=	atoi(argv[5]); // contrast function choice
		unsigned MinParticlesInLeaf = 4; // minimum particles in each leaf of KD Tree
		int TOL_POW                 = atoi(argv[6]); // tolerance for compressions in powers of 10
		int nLevelsUniform = atoi(argv[7]);
		int TOL_POW_NCA_precond = atoi(argv[8]);
		#endif

	#ifdef noToLS
		int sqrtRootN =	atoi(argv[1]);
		double L						=	atof(argv[2]); //domain size
		int TOL_POW         = atoi(argv[3]); // tolerance for compressions in powers of 10
		int nParticlesInLeafAlong1D	=	atoi(argv[4]); // assuming the particles are located at tensor product chebyNodes/uniform
		int TOL_POW_NCA_precond = atoi(argv[5]);
		int nLevelsUniform		=	ceil(log(double(sqrtRootN)/nParticlesInLeafAlong1D)/log(2));
		std::cout << "nLevelsUniform: " << nLevelsUniform << std::endl;
		#endif

	std::vector<pts2D> particles_X, particles_Y; //dummy variables
	#ifdef yestoLS
		userkernel* mykernel		=	new userkernel(nChebNodes, treeAdaptivity, degreeOfBases, L, nLevelsUniform);
		unsigned N = mykernel->N;
		#endif
	#ifdef noToLS
		userkernel* mykernel		=	new userkernel(sqrtRootN, L, nLevelsUniform);
		unsigned N = sqrtRootN*sqrtRootN;
	#endif
	unsigned nLevels            = mykernel->nLevelsResolve;//for quad tree

	double start, end;
	unsigned n_Dimension    =       2;  //      Dimension.
	unsigned n_Properties   =       1;  //      Number of properties/rhs.
	double* locations       =       new double [N*n_Dimension];   //      Stores all the locations.

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
	Vec phi;
	double err;
	//////////// ONES AT RANDOM PLACES  //////////////////////////
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
	#ifdef AIFMM_Direct
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

	aifmm->backSubstitute1(true_Ax); //CONVERGENCE

	start	=	omp_get_wtime();
	aifmm->backSubstitute2(); //CONVERGENCE
	end		=	omp_get_wtime();
	timeFactorize +=	(end-start);

	start	=	omp_get_wtime();
	aifmm->backSubstitute3(); //CONVERGENCE
	end		=	omp_get_wtime();
	double timeSolve =	(end-start);

	aifmm->getPhi(phi);

	err = (true_x - phi).norm()/true_x.norm();

	std::cout << N << " " << aifmm->A->getMaxRank() << " " << timeAssemble << " " << timeFactorize << " " << timeSolve << " " << err << std::endl;
	// exit(0);

	delete aifmm;
	#endif

	//////////////////////////////////////////////////////////////
	#ifdef GMRES
	int TOL_POW_NCA_AFMM = TOL_POW;
	double gmres_tolerance = 1.0e-10;

	start		=	omp_get_wtime();
	AFMM<userkernel>* afmm1 = new AFMM<userkernel>(nParticlesInLeafAlong1D, L, TOL_POW_NCA_AFMM, mykernel->gridPoints, mykernel);
	end		=	omp_get_wtime();
	double timeGMRESInitialise = end-start;
	int maxRank = afmm1->A->getMaxRank();

	classGMRES* G = new classGMRES();
	int maxIterations = 400;
	double errGMRES;
	int noOfIterations;

	std::vector<double> e_gmres;
	start		=	omp_get_wtime();
	G->gmres(afmm1, true_Ax, maxIterations, gmres_tolerance, phi, errGMRES, noOfIterations, e_gmres);
	end		=	omp_get_wtime();
	double timeGMRES = end-start;
	err = (true_x - phi).norm()/true_x.norm();
	std::cout << N << " " << maxRank << " " << noOfIterations << " " << timeGMRESInitialise << " " << timeGMRES << " " << err << std::endl;
	delete afmm1;
	delete G;

	// exit(0);

	#endif

	//// PRECONDITIONER
	#ifdef PRECOND_AIFMM
	int TOL_POW_NCA_AFMM = TOL_POW;
	double gmres_tolerance = 1.0e-12;

	start	=	omp_get_wtime();
	AIFMM<userkernel> *aifmm = new AIFMM<userkernel>(mykernel, N, nLevelsUniform, TOL_POW_NCA_precond, locations);
	end		=	omp_get_wtime();
	double timeAssemble =	(end-start);
	std::cout << std::endl << "Time taken to assemble PRECOND is: " << timeAssemble << std::endl;

	//////////////////////// ELIMINATE PHASE ////////////////////////////////////////
	start	=	omp_get_wtime();
	aifmm->factorize();
	end		=	omp_get_wtime();
	double timeFactorize =	(end-start);
	std::cout << std::endl << "Time taken to factorize PRECOND (eliminate phase): " << timeFactorize << std::endl;
	/////////////////////
	start		=	omp_get_wtime();
	AFMM<userkernel>* afmm1 = new AFMM<userkernel>(nParticlesInLeafAlong1D, L, TOL_POW_NCA_AFMM, mykernel->gridPoints, mykernel);
	end		=	omp_get_wtime();
	double timeGMRESInitialise = end-start;
	std::cout << std::endl << "Time taken by GMRES (AFMM) Initialisation: " << timeGMRESInitialise << std::endl;

	classGMRES* G = new classGMRES();
	int maxIterations = 400;
	double errGMRES;
	int noOfIterations;

	std::cout << std::endl << "begining Hybrid solver...";
	std::vector<double> e_hybrid;
	start		=	omp_get_wtime();
	G->gmres(afmm1, aifmm, true_Ax, maxIterations, gmres_tolerance, phi, errGMRES, noOfIterations, e_hybrid);
	end		=	omp_get_wtime();
	std::cout << "done" << std::endl;
	double timeHybrid = end-start;
	std::cout << std::endl << "Time taken by Hybrid solver: " << timeHybrid << std::endl;
	std::cout << std::endl << "e_Hybrid" << std::endl;
	for (size_t i = 0; i < noOfIterations; i++) {
		std::cout << e_hybrid[i] << " ";
	}
	err = (true_x - phi).norm()/true_x.norm();
	std::cout << std::endl;
	std::cout << std::endl << "timeHybrid: " << timeHybrid << " errHybrid: " << err << "	resGMRES: " << errGMRES << "	noOfIterations: " << noOfIterations << std::endl;

	delete aifmm;
	delete afmm1;
	delete G;
	#endif

	delete locations;
	delete mykernel;
}
