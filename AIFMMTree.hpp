//
//  FMM2DTree.hpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#ifndef _AIFMMTree_HPP__
#define _AIFMMTree_HPP__
#include <fstream>
#include <iostream>
#include <string>
#include <filesystem>
// #include <cereal/archives/portable_binary.hpp>
// #include <cereal/types/vector.hpp>

// namespace storedata
// {
// Eigen::VectorXd inline load_vec(std::string fname)
// {
//         Eigen::VectorXd X;
//         std::vector<double> K;
//         std::ifstream infile (fname, std::ios::binary);
//         cereal::PortableBinaryInputArchive iarchive(infile);      // Create an input archive
//         iarchive(K);  // Read the data from the archive
//         X = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(K.data(), K.size());
//         return X;
// }
//
// void inline save_vec(std::string fname,Eigen::VectorXd X)
// {
//         std::vector<double> K(X.data(), X.data() + X.size());
//         std::ofstream outfile (fname,std::ios::binary);
//         cereal::PortableBinaryOutputArchive oarchive(outfile);    // Create an output archive
//         oarchive(CEREAL_NVP(K));                                  // Write the data to the archive
// }
// }

class FMM2DBox {
public:
	int boxNumber;
	int parentNumber;
	int childrenNumbers[4];
	int neighborNumbers[8];
	int innerNumbers[16];
	int outerNumbers[24];
	bool Eliminated;

	FMM2DBox () {
		boxNumber		=	-1;
		parentNumber	=	-1;
		for (int l=0; l<4; ++l) {
			childrenNumbers[l]	=	-1;
		}
		for (int l=0; l<8; ++l) {
			neighborNumbers[l]	=	-1;
		}
		for (int l=0; l<16; ++l) {
			innerNumbers[l]		=	-1;
		}
		for (int l=0; l<24; ++l) {
			outerNumbers[l]		=	-1;
		}
		Eliminated = false;
	}

	int NumParticles;
	int NumMultipoles; // = NumLocals
	int NumLocals; // = NumLocals
	// int NumLocals;
	pts2D center;
	// Vec charges;
	Vec particles;
	Vec multipoles;
	Vec locals;
	Vec particle_rhs;
	Vec multipole_rhs; //intially value will be 0
	Vec local_rhs; //intially value will be 0

	std::map<int, Mat> M2L;
	std::map<int, Mat> P2P;
	std::map<int, Mat> M2P;
	std::map<int, Mat> P2L;
  std::map<int, Mat> L2P_f;
	// Mat L2P_f;//fill-ins that occur in P2P equation due to elimination of x corresponding to z
	Mat L2M_f; //fill-ins that occur in U.transpose() equation due to elimination of x
	std::map<int, Mat> P2M_f; //fill-ins that occur in U.transpose() equation due to elimination of x
	std::map<int, Mat> M2M_f;
	Mat L2P;					//	Transfer from locals of parent to locals of children.
	Mat P2M;					//	Transfer from multipoles of 4 children to multipoles of parent.
	//	The following will be stored only at the leaf nodes
	Eigen::ColPivHouseholderQR<Mat> U_qr;
	Eigen::ColPivHouseholderQR<Mat> V_qr;
  std::vector<pts2D> particle_loc;
  std::vector<int> chargeLocations;

	std::vector<int> incoming_chargePoints;//equivalent points {y_{k}^{B,i}}
	std::vector<int> incoming_checkPoints;//check points {x_{k}^{B,i}}
	std::vector<int> outgoing_chargePoints;
	std::vector<int> outgoing_checkPoints;
	Mat Ac;
	Mat Ar;
};

template <typename kerneltype>
class FMM2DTree {
public:
	kerneltype* K;
	int nLevels;			//	Number of levels in the tree.
	int N;					//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	double smallestBoxSize;	//	This is L/2.0^(nLevels).
	// double a;				//	Cut-off for self-interaction. This is less than the length of the smallest box size.

	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> boxRadius;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<std::vector<FMM2DBox> > tree;	//	The tree storing all the information.

	int nParticlesInLeafAlong1D;
  int nParticlesInLeaf;
  //	Chebyshev nodes
	std::vector<double> Nodes1D;
	std::vector<pts2D> Nodes;
  std::vector<pts2D> gridPoints; //all particles in domain
	int TOL_POW;
	double RRQR_threshold,LU_threshold;
	Vec b_error_check;
	double* locations;
	// std::vector<int> boxNumbers;
	// std::vector<int> NumberOfParticlesInLeaves;
	std::vector<std::pair<int, int> > P2L_M2P;
	std::vector<std::pair<int, int> > P2P;
// public:
	FMM2DTree(kerneltype* K, int N, int nLevels, int TOL_POW, double* locations) {
		this->K					=	K;
		this->nLevels		=	nLevels;
		this->L					=	1.0;
		this->locations = locations;
		// this->boxNumbers= boxNumbers;
		// this->NumberOfParticlesInLeaves = NumberOfParticlesInLeaves;
    // this->nParticlesInLeafAlong1D = nParticlesInLeafAlong1D;
    // this->nParticlesInLeaf = nParticlesInLeafAlong1D*nParticlesInLeafAlong1D;
		this->TOL_POW = TOL_POW;
		this->RRQR_threshold = pow(10,-1.0*TOL_POW);
		this->LU_threshold = pow(10,-1.0*TOL_POW);
    nBoxesPerLevel.push_back(1);
		boxRadius.push_back(L);
		for (int k=1; k<=nLevels; ++k) {
			nBoxesPerLevel.push_back(4*nBoxesPerLevel[k-1]);
			boxRadius.push_back(0.5*boxRadius[k-1]);
		}
		this->smallestBoxSize	=	boxRadius[nLevels];
		K->a					=	smallestBoxSize;
		this->N				=	N;
	}

	void createTree() {
		//	First create root and add to tree
		FMM2DBox root;
		root.boxNumber		=	0;
		root.parentNumber	=	-1;
		#pragma omp parallel for
		for (int l=0; l<4; ++l) {
			root.childrenNumbers[l]	=	l;
		}
		#pragma omp parallel for
		for (int l=0; l<8; ++l) {
			root.neighborNumbers[l]	=	-1;
		}
		#pragma omp parallel for
		for (int l=0; l<16; ++l) {
			root.innerNumbers[l]	=	-1;
		}
		#pragma omp parallel for
		for (int l=0; l<24; ++l) {
			root.outerNumbers[l]	=	-1;
		}
		std::vector<FMM2DBox> rootLevel;
		rootLevel.push_back(root);
		tree.push_back(rootLevel);

		for (int j=1; j<=nLevels; ++j) {
			std::vector<FMM2DBox> level;
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				FMM2DBox box;
				box.boxNumber		=	k;
				box.parentNumber	=	k/4;
				for (int l=0; l<4; ++l) {
					box.childrenNumbers[l]	=	4*k+l;
				}
				level.push_back(box);
			}
			tree.push_back(level);
		}
	}

	//	Assigns the interactions for child0 of a box
	void assign_Child0_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N5  |  N4  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  **  |  N3  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[3]	=	nC+1;
			tree[nL][nC].neighborNumbers[4]	=	nC+2;
			tree[nL][nC].neighborNumbers[5]	=	nC+3;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/****************************/
			/*				   ______	*/
			/*				  |		 |	*/
			/*				  |	 **  |	*/
			/*	 _____________|______|  */
			/*	|	   |	  |			*/
			/*	|  I15 |  N0  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I0  |  I1  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 ______			  		*/
			/*	|	   |	  			*/
			/*	|  **  |				*/
			/*	|______|______			*/
			/*	|	   |	  |			*/
			/*	|  N1  |  N2  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I2  |  I3  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's second neighbor
		{
			/************************************/
			/*	 ______			  				*/
			/*	|	   |	  					*/
			/*	|  **  |						*/
			/*	|______|	   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  I5  |  O8  |	*/
			/*				  |______|______|	*/
			/*				  |	     |	    |	*/
			/*				  |  I4  |  O7  |	*/
			/*				  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[8]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's third neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  I7  |  O10 |	*/
			/*	 ______		  |______|______|	*/
			/*	|	   |	  |	     |	    |	*/
			/*	|  **  |	  |  I6  |  O9  |	*/
			/*	|______|	  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN!=-1) {
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[9]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[3];
			}

		}

		//	Assign children of parent's fourth neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  O13 |  O12 |	*/
			/*				  |______|______|	*/
			/*				  |	     |	    |	*/
			/*		    	  |  I8  |  O11 |	*/
			/*		    	  |______|______|	*/
			/*									*/
			/*									*/
			/*	 ______							*/
			/*  |      |						*/
			/*  |  **  |						*/
			/*  |______|						*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[11]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[12]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[13]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  O15 |  O14 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  I10 |  I9  |		*/
			/*	|______|______|		*/
			/*						*/
			/*						*/
			/*	 ______				*/
			/*  |	   |			*/
			/*	|  **  |			*/
			/*	|______|			*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[14]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[15]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/****************************/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  O17 |  O16 |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I12 |  I11 |			*/
			/*	|______|______|			*/
			/*							*/
			/*							*/
			/*				   ______	*/
			/*  			  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[16]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[17]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/****************************/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I13 |  N6  |			*/
			/*	|______|______|______	*/
			/*  |	   |	  |		 |	*/
			/*	|  I14 |  N7  |	 **  |	*/
			/*	|______|______|______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[3];
			}
		}
	}

	//	Assigns the interactions for child1 of a box
	void assign_Child1_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+1;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N6  |  N5  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N7  |  **  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[7]	=	nC-1;
			tree[nL][nC].neighborNumbers[5]	=	nC+1;
			tree[nL][nC].neighborNumbers[6]	=	nC+2;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/************************************/
			/*				   		  ______	*/
			/*				  	     |		|	*/
			/*				         |	**  |	*/
			/*	 _____________       |______|  	*/
			/*	|	   |	  |					*/
			/*	|  O22 |  I15 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O23 |  I0  |					*/
			/*	|______|______|					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[23]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[22]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 		______		  	*/
			/*		   |	  |			*/
			/*	       |  **  |			*/
			/*	 ______|______|			*/
			/*	|	   |	  |			*/
			/*	|  N0  |  N1  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I1  |  I2  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's second neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |	  			*/
			/*	|______|_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  N2  |  I5  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*	       |  I3  |  I4  |	*/
			/*		   |______|______|	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's third neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  N4  |	 I7	 |  */
			/*	 ______|______|______|	*/
			/*	|	   |	  |	     |	*/
			/*	|  **  |  N3  |  I6  |  */
			/*	|______|______|______| 	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN != -1){
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  O14 |  O13 |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  I9  |  I8  |	*/
			/*		   |______|______|	*/
			/*				  			*/
			/*				  			*/
			/*	 ______					*/
			/*	|	   |				*/
			/*  |  **  |				*/
			/*  |______|				*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[13]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[14]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  O16 |  O15 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  I11 |  I10 |		*/
			/*	|______|______|		*/
			/*						*/
			/*						*/
			/*		    ______		*/
			/* 		   |	  |		*/
			/*		   |  **  |		*/
			/*		   |______|		*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[15]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[16]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/************************************/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O18 |  O17 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O19 |  I12 |					*/
			/*	|______|______|					*/
			/*									*/
			/*									*/
			/*				   		  ______	*/
			/*  			  		 |		|	*/
			/*				  		 |	** 	|	*/
			/*				  		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[19]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[17]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[18]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/************************************/
			/*									*/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O20 |  I13 |					*/
			/*	|______|______|		  ______	*/
			/*  |	   |	  |		 |		|	*/
			/*	|  O21 |  I14 |	 	 |	**  |	*/
			/*	|______|______|		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[21]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[20]	=	tree[j][nN].childrenNumbers[3];
			}
		}
	}

	//	Assigns the interactions for child2 of a box
	void assign_Child2_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+2;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N7  |  **  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N0  |  N1  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[0]	=	nC-2;
			tree[nL][nC].neighborNumbers[1]	=	nC-1;
			tree[nL][nC].neighborNumbers[7]	=	nC+1;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/************************************/
			/*				   		  ______	*/
			/*				  	     |		|	*/
			/*				         |	**  |	*/
			/*				         |______|  	*/
			/*									*/
			/*									*/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O23 |  I0  |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O0  |  O1  |					*/
			/*	|______|______|					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[23]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 		______		  	*/
			/*		   |	  |			*/
			/*	       |  **  |			*/
			/*	 	   |______|			*/
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I1  |  I2  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O2  |  O3  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's second neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |	  			*/
			/*	|______|				*/
			/*							*/
			/*							*/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  I3  |  I4  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*	       |  O4  |  O5  |	*/
			/*		   |______|______|	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[5]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's third neighbor
		{
			/****************************/
			/*	 ____________________	*/
			/*	|	   |	  |	     |	*/
			/*	|  **  |  N3  |	 I6	 |  */
			/*	|______|______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  N2  |  I5  |  */
			/*		   |______|______| 	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  I9  |  I8  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  N4  |  I7  |	*/
			/*	 ______|______|______|	*/
			/*	|	   |	  			*/
			/*	|  **  |	  			*/
			/*	|______|	  			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  I11 |  I10 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N6  |  N5  |		*/
			/*	|______|______|		*/
			/*		   |	  |		*/
			/*		   |  **  |		*/
			/*		   |______|		*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/************************************/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O19 |  I12 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O20 |  I13 |					*/
			/*	|______|______|		  ______	*/
			/*  			  		 |		|	*/
			/*				  		 |	** 	|	*/
			/*				  		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[20]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[19]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/************************************/
			/*									*/
			/*	 _____________		  ______	*/
			/*	|	   |	  |		 |	    |	*/
			/*	|  O21 |  I14 |		 |	**	|	*/
			/*	|______|______|		 |______|	*/
			/*  |	   |	  |		 			*/
			/*	|  O22 |  I15 |	 	 			*/
			/*	|______|______|		 			*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[22]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[21]	=	tree[j][nN].childrenNumbers[3];
			}
		}
	}

	//	Assigns the interactions for child3 of a box
	void assign_Child3_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+3;
		int nN;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  **  |  N3  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N1  |  N2  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[1]	=	nC-3;
			tree[nL][nC].neighborNumbers[2]	=	nC-2;
			tree[nL][nC].neighborNumbers[3]	=	nC-1;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/****************************/
			/*				   ______	*/
			/*				  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|  */
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I0  |  I1  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O1  |  O2  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[0];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |				*/
			/*	|______|				*/
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I2  |  I3  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O3  |  O4  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's second neighbor
		{
			/************************************/
			/*	 ______		  					*/
			/*	|	   |						*/
			/*	|  **  |	  					*/
			/*	|______|						*/
			/*									*/
			/*									*/
			/*				   _____________	*/
			/*		   		  |	     |	    |	*/
			/*		   		  |  I4  |  O7  |	*/
			/*		   		  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*	       		  |  O5  |  O6  |	*/
			/*		   		  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[2];
			if (nN != -1) {
				tree[nL][nC].outerNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's third neighbor
		{
			/************************************/
			/*	 ______		   _____________	*/
			/*	|	   |	  |	     |		|	*/
			/*	|  **  |      |	 I6	 |  O9	|	*/
			/*	|______|	  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I5  |  O8  |  	*/
			/*		   		  |______|______| 	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[8]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[9]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I8  |  O11 |	*/
			/*		   		  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I7  |  O10 |	*/
			/*	 ______	      |______|______|	*/
			/*	|	   |	  					*/
			/*	|  **  |	  					*/
			/*	|______|	  					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[11]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  I10 |  I9  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N5  |  N4  |		*/
			/*	|______|______|		*/
			/*	|	   |			*/
			/*	|  **  |			*/
			/*	|______|			*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/****************************/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I12 |  I11 |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I13 |  N6  |			*/
			/*	|______|______|______	*/
			/*  			  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[6];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[3];
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/****************************/
			/*							*/
			/*	 ____________________	*/
			/*	|	   |	  |		 |	*/
			/*	|  I14 |  N7  |	 **	 |	*/
			/*	|______|______|______|	*/
			/*  |	   |	  |		 	*/
			/*	|  I15 |  N0  |	 	 	*/
			/*	|______|______|		 	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[3];
			}
		}
	}


	//	Assigns the interactions for the children of a box
	void assign_Box_Interactions(int j, int k) {
		assign_Child0_Interaction(j,k);
		assign_Child1_Interaction(j,k);
		assign_Child2_Interaction(j,k);
		assign_Child3_Interaction(j,k);
	}

	//	Assigns the interactions for the children all boxes at a given level
	void assign_Level_Interactions(int j) {
		#pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			assign_Box_Interactions(j,k);
		}
	}

	//	Assigns the interactions for the children all boxes in the tree
	void assign_Tree_Interactions() {
		for (int j=0; j<nLevels; ++j) {
			assign_Level_Interactions(j);
		}
	}

	void assign_Center_Location() {
		int J;
		tree[0][0].center.x	=	0.0;
		tree[0][0].center.y	=	0.0;
		for (int j=0; j<nLevels; ++j) {
			J	=	j+1;
			double shift	=	0.5*boxRadius[j];
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				tree[J][4*k].center.x		=	tree[j][k].center.x-shift;
				tree[J][4*k+1].center.x	=	tree[j][k].center.x+shift;
				tree[J][4*k+2].center.x	=	tree[j][k].center.x+shift;
				tree[J][4*k+3].center.x	=	tree[j][k].center.x-shift;

				tree[J][4*k].center.y		=	tree[j][k].center.y-shift;
				tree[J][4*k+1].center.y	=	tree[j][k].center.y-shift;
				tree[J][4*k+2].center.y	=	tree[j][k].center.y+shift;
				tree[J][4*k+3].center.y	=	tree[j][k].center.y+shift;
			}
		}
	}

	// void assignLeafChargeLocations() {
	// 	for (size_t i = 0; i < N*2; i+=2) {
	// 		pts2D temp;
	// 		temp.x = locations[i];
	// 		temp.y = locations[i+1];
	// 		gridPoints.push_back(temp);
	// 	}
	// 	int startIndex = 0;
	// 	for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
	// 		int boxNum = boxNumbers[k];
	// 		int NumParticles = NumberOfParticlesInLeaves[k];
	// 		for (size_t i = 0; i < NumParticles; i++) {
	// 			tree[nLevels][boxNum].chargeLocations.push_back(startIndex+i);
	// 		}
	// 		startIndex += NumParticles;
	// 	}
	// 	K->particles_X = gridPoints;//object of base class FMM_Matrix
	// 	K->particles_Y = gridPoints;
	// }
/////
void reorder(Vec &potential) {
	Vec potentialTemp = potential;
	int start = 0;
	for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
		for (size_t i = 0; i < tree[nLevels][k].chargeLocations.size(); i++) {
			int index = tree[nLevels][k].chargeLocations[i];
			potential(index) = potentialTemp(start);
			start++;
		}
	}
}

	void assignChargeLocations() {
		// K->particles = Nodes;//object of base class FMM_Matrix
		for (size_t i = 0; i < N*2; i+=2) {
			pts2D temp;
			temp.x = locations[i];
			temp.y = locations[i+1];
			gridPoints.push_back(temp);
		}
		// std::cout << "N: " << N << std::endl;
		// for (size_t i = 0; i < N; i++) {
		// 	std::cout << gridPoints[i].x << " " << gridPoints[i].y << std::endl;
		// }
		K->particles_X = gridPoints;//object of base class FMM_Matrix
		K->particles_Y = gridPoints;

		for (size_t i = 0; i < N; i++) {
			tree[0][0].chargeLocations.push_back(i);
		}
		for (size_t j = 0; j < nLevels; j++) { //assign particles to its children
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				int J = j+1;
				int Kp = 4*k;
				for (size_t i = 0; i < tree[j][k].chargeLocations.size(); i++) {
					int index = tree[j][k].chargeLocations[i];
					if (gridPoints[index].x <= tree[j][k].center.x) { //children 0,3
						if (gridPoints[index].y <= tree[j][k].center.y) { //child 0
							tree[J][Kp].chargeLocations.push_back(index);
						}
						else { //child 3
							tree[J][Kp+3].chargeLocations.push_back(index);
						}
					}
					else { //children 1,2
						if (gridPoints[index].y <= tree[j][k].center.y) { //child 1
							tree[J][Kp+1].chargeLocations.push_back(index);
						}
						else { //child 2
							tree[J][Kp+2].chargeLocations.push_back(index);
						}
					}
				}
			}
		}
	}

	int getMaxRank() {
		int max = 0;
		for (int j = nLevels; j >= 2; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (max < tree[j][k].multipoles.size()) {
					max = tree[j][k].multipoles.size();
				}
			}
		}
		return max;
	}

	void check55() {
		for (int j = nLevels; j >= 2; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				std::cout << "j: " << j << "	k: " << k << "	cL: " << tree[j][k].chargeLocations.size() << std::endl;
				for (size_t i = 0; i < tree[j][k].chargeLocations.size(); i++) {
					std::cout << tree[j][k].chargeLocations[i] << ",";
				}
				std::cout << std::endl;
				// tree[j][k].chargeLocations.clear();
			}
		}
		for (int j = nLevels; j >= 2; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				std::cout << "j: " << j << "	k: " << k << "	pr: " << tree[j][k].particle_rhs.size() << std::endl;
				for (size_t i = 0; i < tree[j][k].particle_rhs.size(); i++) {
					std::cout << tree[j][k].particle_rhs[i] << ",";
				}
				std::cout << std::endl;
				// tree[j][k].chargeLocations.clear();
			}
		}
	}

	void assignNonLeafChargeLocations() {
		for (int j = nLevels-1; j >= 1; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				tree[j][k].chargeLocations.clear();
				// for (size_t c = 0; c < 4; c++) {
				// 	tree[j][k].chargeLocations.insert(tree[j][k].chargeLocations.end(), tree[j+1][4*k+c].chargeLocations.begin(), tree[j+1][4*k+c].chargeLocations.end());
				// }
				// std::cout << "tree[j][k].charges: " << tree[j][k].charges.size() << std::endl;
			}
		}
	}
//////

	// void assign_Leaf_rhs(Vec &rhs) {
	// 	int index = 0;
	// 	for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
	// 		int boxNum = boxNumbers[k];
	// 		int NumParticles = NumberOfParticlesInLeaves[k];
	// 		tree[nLevels][boxNum].particle_rhs = Vec(NumParticles);
	// 		for (size_t i = 0; i < NumParticles; i++) {
	// 			tree[nLevels][boxNum].particle_rhs(i) = rhs(index+i);
	// 		}
	// 		index += NumParticles;
	// 	}
	// }

	void assign_Leaf_rhs(Vec &charges) {
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			tree[nLevels][k].particle_rhs	=	Vec::Zero(tree[nLevels][k].chargeLocations.size());
			for (size_t i = 0; i < tree[nLevels][k].chargeLocations.size(); i++) {
				int index = tree[nLevels][k].chargeLocations[i];
				tree[nLevels][k].particle_rhs[i]	=	charges[index];
			}
		}
	}


	void get_L2P_P2M_box(int j, int k) {
		if (tree[j][k].incoming_checkPoints.size() < tree[j][k].outgoing_checkPoints.size()) {
			tree[j][k].outgoing_checkPoints.erase(tree[j][k].outgoing_checkPoints.begin()+tree[j][k].incoming_checkPoints.size(), tree[j][k].outgoing_checkPoints.end());
			tree[j][k].outgoing_chargePoints.erase(tree[j][k].outgoing_chargePoints.begin()+tree[j][k].incoming_checkPoints.size(), tree[j][k].outgoing_chargePoints.end());
		}
		else {
			tree[j][k].incoming_checkPoints.erase(tree[j][k].incoming_checkPoints.begin()+tree[j][k].outgoing_checkPoints.size(), tree[j][k].incoming_checkPoints.end());
			tree[j][k].incoming_chargePoints.erase(tree[j][k].incoming_chargePoints.begin()+tree[j][k].outgoing_checkPoints.size(), tree[j][k].incoming_chargePoints.end());
		}
		// std::cout << "here " << std::endl;
		Mat temp1 = tree[j][k].Ac.block(0,0,tree[j][k].Ac.rows(), tree[j][k].incoming_checkPoints.size());
		tree[j][k].L2P = Mat(tree[j][k].Ac.rows(),tree[j][k].incoming_checkPoints.size());
		// std::cout << "here " << std::endl;
		if (tree[j][k].incoming_checkPoints.size() > 0) {
			Mat D1 = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][k].incoming_chargePoints);
			Eigen::PartialPivLU<Mat> D1_lu = Eigen::PartialPivLU<Mat>(D1.adjoint());
			tree[j][k].L2P = D1_lu.solve(temp1.adjoint()).adjoint();
		}
		// std::cout << "here " << std::endl;
		Mat temp2 = tree[j][k].Ar.block(0,0,tree[j][k].outgoing_checkPoints.size(),tree[j][k].Ar.cols());
		tree[j][k].P2M = Mat(tree[j][k].outgoing_checkPoints.size(),tree[j][k].Ar.cols());
		// std::cout << "here " << std::endl;
		if (tree[j][k].outgoing_checkPoints.size() > 0) {
			Mat D2 = K->getMatrix(tree[j][k].outgoing_checkPoints, tree[j][k].outgoing_chargePoints);
			Eigen::PartialPivLU<Mat> D2_lu = Eigen::PartialPivLU<Mat>(D2);
			tree[j][k].P2M = D2_lu.solve(temp2);
		}

		// std::cout << "j: " << j << "	k: " << k << "	N: " << abs(tree[j][k].L2P.rows()-tree[j][k].P2M.cols()) << "	p: " << abs(tree[j][k].P2M.rows() - tree[j][k].L2P.cols()) << std::endl;
		// std::cout << "j: " << j << "	k: " << k << "	I: " << tree[j][k].incoming_checkPoints.size() << "	O: " << tree[j][k].outgoing_checkPoints.size() << "	L2P.N: " << tree[j][k].L2P.norm()  << "	P2M.N: " << tree[j][k].P2M.norm() << std::endl;// << ", " << tree[j][k].L2P.cols() << "	P2M: " << tree[j][k].P2M.rows() << ", " << tree[j][k].P2M.cols() << std::endl;
	}

	void check1() {
		for (int j=nLevels; j>=2; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				std::cout << "j: " << j << "	k: " << k << "	L2P: " << tree[j][k].L2P.rows() << ", " << tree[j][k].L2P.cols() << "	P2M: " << tree[j][k].P2M.rows() << ", " << tree[j][k].P2M.cols() << "	D: " << tree[j][k].L2P.rows() - tree[j][k].P2M.cols() << ", " << tree[j][k].L2P.cols() - tree[j][k].P2M.rows() << std::endl;
			}
		}
	}

	void get_L2P_P2M_level(int j) {
		#pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			get_L2P_P2M_box(j, k);
		}
	}

	void getNodes() {
		for (int j=nLevels; j>=2; j--) {
			getNodes_outgoing_level(j);
			getNodes_incoming_level(j);
			get_L2P_P2M_level(j);
		}
	}

	void getNodes_outgoing_level(int j) { //LFR; box interactions
		// int ComputedRank, n_rows, n_cols;
    #pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			getNodes_outgoing_box(j, k);
      // getNodes_outgoing_box(j, k, n_rows, n_cols, ComputedRank);
		}
	}

	void getNodes_incoming_level(int j) { //LFR; box interactions
		// int ComputedRank, n_rows, n_cols;
    #pragma omp parallel for
    for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			getNodes_incoming_box(j, k);
      // getNodes_incoming_box(j, k, n_rows, n_cols, ComputedRank);
		}
	}

	void getParticlesFromChildren_incoming_row(int j, int k, std::vector<int>& searchNodes) {
		if (j==nLevels) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			for (int c = 0; c < 4; c++) {
				searchNodes.insert(searchNodes.end(), tree[J][4*k+c].incoming_checkPoints.begin(), tree[J][4*k+c].incoming_checkPoints.end());
			}
		}
	}

	void getParticlesFromChildren_incoming_col(int j, int k, std::vector<int>& searchNodes) {
		if (j==nLevels) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			for (int c = 0; c < 4; c++) {
				searchNodes.insert(searchNodes.end(), tree[J][4*k+c].outgoing_chargePoints.begin(), tree[J][4*k+c].outgoing_chargePoints.end());
			}
		}
	}

	void getParticlesFromChildren_outgoing_row(int j, int k, std::vector<int>& searchNodes) {
		if (j==nLevels) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			for (int c = 0; c < 4; c++) {
				searchNodes.insert(searchNodes.end(), tree[J][4*k+c].incoming_checkPoints.begin(), tree[J][4*k+c].incoming_checkPoints.end());
			}
		}
	}

	void getParticlesFromChildren_outgoing_col(int j, int k, std::vector<int>& searchNodes) {
		if (j==nLevels) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			for (int c = 0; c < 4; c++) {
				searchNodes.insert(searchNodes.end(), tree[J][4*k+c].outgoing_chargePoints.begin(), tree[J][4*k+c].outgoing_chargePoints.end());
			}
		}
	}

	void getNodes_outgoing_box(int j, int k) {
  // void getNodes_outgoing_box(int j, int k, int& n_rows, int& n_cols, int& ComputedRank) {
		int n_rows, n_cols, ComputedRank;
    std::vector<int> boxA_Nodes;
		getParticlesFromChildren_outgoing_col(j, k, boxA_Nodes);
		std::vector<int> IL_Nodes;//indices
		for(int in=0; in<16; ++in) {
			if(tree[j][k].innerNumbers[in] != -1) {
				int kIL = tree[j][k].innerNumbers[in];
				std::vector<int> chargeLocations;
				getParticlesFromChildren_outgoing_row(j, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}
		}
		for(int on=0; on<24; ++on) {
			if(tree[j][k].outerNumbers[on] != -1) {
				int kIL = tree[j][k].outerNumbers[on];
				std::vector<int> chargeLocations;
				getParticlesFromChildren_outgoing_row(j, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}
		}
		n_rows = IL_Nodes.size();
		n_cols = boxA_Nodes.size();
		// int tol_pow = TOL_POW;
		// double tol_ACA = pow(10,-1.0*tol_pow);
		std::vector<int> row_indices, col_indices;
		row_indices = IL_Nodes;
		col_indices = boxA_Nodes;//object of base class G_LowRank
		std::vector<int> row_bases, col_bases;
		Mat Ac;
		LowRank* LR		=	new LowRank(K, pow(10,-TOL_POW), row_indices, col_indices);
		LR->rookPiv(row_bases, col_bases, ComputedRank, Ac, tree[j][k].Ar);
		// LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, tree[j][k].Ar);

		int minN = n_rows;
		if (n_rows > n_cols) {
			minN = n_cols;
		}
		if(ComputedRank > 0) {
			for (int r = 0; r < row_bases.size(); r++) {
				tree[j][k].outgoing_checkPoints.push_back(IL_Nodes[row_bases[r]]);
			}
			for (int c = 0; c < col_bases.size(); c++) {
				tree[j][k].outgoing_chargePoints.push_back(boxA_Nodes[col_bases[c]]);
			}
			// Mat D = K->getMatrix(tree[j][k].outgoing_checkPoints, tree[j][k].outgoing_chargePoints);
			// Eigen::PartialPivLU<Mat> D_lu = Eigen::PartialPivLU<Mat>(D);
			// tree[j][k].P2M = D_lu.solve(tree[j][k].Ar);
		}
		// else {
		// 	tree[j][k].P2M = Mat(tree[j][k].outgoing_checkPoints.size(),tree[j][k].Ar.cols());
		// }
		// std::cout << "j: " << j << "	k: "<< k << "	ComputedRank: " << ComputedRank << std::endl;
	}

	void getNodes_incoming_box(int j, int k) {
  // void getNodes_incoming_box(int j, int k, int& n_rows, int& n_cols, int& ComputedRank) {
		int n_rows, n_cols, ComputedRank;
    std::vector<int> boxA_Nodes;
		getParticlesFromChildren_incoming_row(j, k, boxA_Nodes);
		std::vector<int> IL_Nodes;//indices
		for(int in=0; in<16; ++in) {
			if(tree[j][k].innerNumbers[in] != -1) {
				int kIL = tree[j][k].innerNumbers[in];
				std::vector<int> chargeLocations;
				getParticlesFromChildren_incoming_col(j, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}
		}
		for(int on=0; on<24; ++on) {
			if(tree[j][k].outerNumbers[on] != -1) {
				int kIL = tree[j][k].outerNumbers[on];
				std::vector<int> chargeLocations;
				getParticlesFromChildren_incoming_col(j, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}
		}
		n_rows = boxA_Nodes.size();
		n_cols = IL_Nodes.size();
		// int tol_pow = TOL_POW;
		// double tol_ACA = pow(10,-1.0*tol_pow);
		std::vector<int> row_indices, col_indices;
		row_indices = boxA_Nodes;//object of base class G_LowRank
		col_indices = IL_Nodes;
		std::vector<int> row_bases, col_bases;
		Mat Ar;
		LowRank* LR		=	new LowRank(K, pow(10,-TOL_POW), row_indices, col_indices);
		// LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, tree[j][k].Ac, Ar);
		LR->rookPiv(row_bases, col_bases, ComputedRank, tree[j][k].Ac, Ar);

		// std::cout << "j: " << j << "	k: " << k << "	I ComputedRank: " << ComputedRank << std::endl;
		// std::cout << "row_bases: " << row_bases.size() << "	col_bases: " << col_bases.size() << "	tree[j][k].Ac: " << tree[j][k].Ac.rows() << "," << tree[j][k].Ac.cols() << "	Ar: " << Ar.rows() << "," << Ar.cols() << std::endl;
		int minN = n_rows;
		if (n_rows > n_cols) {
			minN = n_cols;
		}
		// tree[j][k].L2P = tree[j][k].Ac;
		if(ComputedRank > 0) {
			for (int r = 0; r < row_bases.size(); r++) {
				tree[j][k].incoming_checkPoints.push_back(boxA_Nodes[row_bases[r]]);
			}
			for (int c = 0; c < col_bases.size(); c++) {
				tree[j][k].incoming_chargePoints.push_back(IL_Nodes[col_bases[c]]);
			}
			// Mat D = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][k].incoming_chargePoints);
			// Eigen::PartialPivLU<Mat> D_T = Eigen::PartialPivLU<Mat>(D.transpose());
			// tree[j][k].L2P = D_T.solve(tree[j][k].Ac.transpose()).transpose();
		}
		// else {
		// 	tree[j][k].L2P = Mat(tree[j][k].Ac.rows(), tree[j][k].incoming_checkPoints.size());
		// }
    // std::cout << "j: " << j << "	k: "<< k << "	ComputedRank: " << ComputedRank << std::endl;
	}

	void assemble_M2L() {
		// #pragma omp parallel for
		for (size_t j = 2; j <= nLevels; j++) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				// #pragma omp parallel for
				for(int in=0; in<16; ++in) {
					if(tree[j][k].innerNumbers[in] != -1) {
						int kIL = tree[j][k].innerNumbers[in];
						// tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].incoming_checkPoints);
						tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].outgoing_chargePoints);
					}
				}
				// #pragma omp parallel for
				for(int on=0; on<24; ++on) {
					if(tree[j][k].outerNumbers[on] != -1) {
						int kIL = tree[j][k].outerNumbers[on];
						// tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].incoming_checkPoints);
						tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].outgoing_chargePoints);
					}
				}
			}
		}
	}

	void check3() {
		std::cout << "is_well_separated(3,6,10): " << is_well_separated(3,6,10) << std::endl;
		std::cout << "is_well_separated(3,7,11): " << is_well_separated(3,7,11) << std::endl;
		std::cout << "tree[3][6].M2L[10]: " << std::endl << tree[3][6].M2L[10] << std::endl << std::endl;
		std::cout << "tree[3][7].M2L[11]: " << std::endl << tree[3][7].M2L[11] << std::endl << std::endl;
	}

	void make_basis_unitary() {
		for (size_t j = 2; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				Eigen::ColPivHouseholderQR<Mat> qr(tree[j][k].L2P.rows(), tree[j][k].L2P.cols());
				tree[j][k].U_qr = qr.compute(tree[j][k].L2P);
				Mat Q_u = qr.householderQ();
				if(tree[j][k].L2P.rows() > tree[j][k].L2P.cols()) {
					Q_u = Q_u.block(0,0,tree[j][k].L2P.rows(), tree[j][k].L2P.cols());
				}
				tree[j][k].L2P = Q_u;
			}
		}
	}

	void update_M2L_to_unitary_basis() {
		// #pragma omp parallel for
		for (size_t j = 2; j <= nLevels; j++) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				// #pragma omp parallel for
				for(int in=0; in<16; ++in) {
					if(tree[j][k].innerNumbers[in] != -1) {
						int kIL = tree[j][k].innerNumbers[in];
						// if (tree[j][k].M2L[kIL].size() == 0) {
							Mat R_u = tree[j][k].U_qr.matrixQR().triangularView<Eigen::Upper>();
							Mat R_v = tree[j][kIL].U_qr.matrixQR().triangularView<Eigen::Upper>();
							if(tree[j][k].L2P.rows() > tree[j][k].L2P.cols()) {
								R_u = R_u.block(0,0,tree[j][k].L2P.cols(), tree[j][k].L2P.cols());
							}
							if(tree[j][kIL].L2P.rows() > tree[j][kIL].L2P.cols()) {
								R_v = R_v.block(0,0,tree[j][kIL].L2P.cols(), tree[j][kIL].L2P.cols());
							}
							tree[j][k].M2L[kIL] = R_u * tree[j][k].M2L[kIL] * R_v.transpose();
						// }
						// if (tree[j][kIL].M2L[k].size() == 0)
						// 	tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose();
					}
				}
				// #pragma omp parallel for
				for(int on=0; on<24; ++on) {
					if(tree[j][k].outerNumbers[on] != -1) {
						int kIL = tree[j][k].outerNumbers[on];
						// if (tree[j][k].M2L[kIL].size() == 0) {
							Mat R_u = tree[j][k].U_qr.matrixQR().triangularView<Eigen::Upper>();
							Mat R_v = tree[j][kIL].U_qr.matrixQR().triangularView<Eigen::Upper>();
							if(tree[j][k].L2P.rows() > tree[j][k].L2P.cols()) {
								R_u = R_u.block(0,0,tree[j][k].L2P.cols(), tree[j][k].L2P.cols());
							}
							if(tree[j][kIL].L2P.rows() > tree[j][kIL].L2P.cols()) {
								R_v = R_v.block(0,0,tree[j][kIL].L2P.cols(), tree[j][kIL].L2P.cols());
							}
							tree[j][k].M2L[kIL] = R_u * tree[j][k].M2L[kIL] * R_v.transpose();
						// }
						// if (tree[j][kIL].M2L[k].size() == 0)
						// 	tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose();
					}
				}
			}
		}
	}

		void compress_P2P_Two_Ways(int j, int k, int b) {
			Mat new_L2P_k;
			Mat new_P2M_transpose_b;
			Mat new_L2P_k_modified;
			Mat new_P2M_transpose_b_modified;
			Mat M2L_k_b;
			compress_P2P(j, k, b, new_L2P_k, new_P2M_transpose_b, M2L_k_b);

			Mat new_L2P_b;
			Mat new_P2M_transpose_k;
			Mat new_L2P_b_modified;
			Mat new_P2M_transpose_k_modified;
			Mat M2L_b_k;
			compress_P2P(j, b, k, new_L2P_b, new_P2M_transpose_k, M2L_b_k);

			if (new_L2P_k.cols() < new_P2M_transpose_k.cols()) {
				new_P2M_transpose_k_modified = new_P2M_transpose_k.block(0,0,new_P2M_transpose_k.rows(),new_L2P_k.cols());
				new_L2P_k_modified = new_L2P_k;
			}
			else {
				new_L2P_k_modified = new_L2P_k.block(0,0,new_L2P_k.rows(),new_P2M_transpose_k.cols());
				new_P2M_transpose_k_modified = new_P2M_transpose_k;
			}
			// if (j==3 && k==19 && b==23) {
			// 	std::cout << "Vt:" << std::endl << new_P2M_transpose_b << std::endl;
			// 	std::cout << "VVt:" << std::endl << new_P2M_transpose_b.adjoint()*new_P2M_transpose_b << std::endl;
			// }

			if (new_L2P_b.cols() < new_P2M_transpose_b.cols()) {
				new_P2M_transpose_b_modified = new_P2M_transpose_b.block(0,0,new_P2M_transpose_b.rows(),new_L2P_b.cols());
				new_L2P_b_modified = new_L2P_b;
				// if (j==3 && k==19 && b==23) {
				// 	std::cout << "Vt: " << std::endl << new_P2M_transpose_b.block(0,0,new_P2M_transpose_b.rows(),new_L2P_b.cols()) << std::endl << std::endl;
				// 	exit(0);
				// }
			}
			else {
				new_L2P_b_modified = new_L2P_b.block(0,0,new_L2P_b.rows(),new_P2M_transpose_b.cols());
				new_P2M_transpose_b_modified = new_P2M_transpose_b;
			}
			tree[j][k].M2L[b] = M2L_k_b.block(0,0,new_L2P_k_modified.cols(),new_P2M_transpose_b_modified.cols());
			// if (j==3 && k==6 && b==10) {
			// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree[j][k].M2L[b]: " << std::endl << tree[j][k].M2L[b] << std::endl << std::endl;
			// }
			update_other_M2Ls_at_k(j, k, b, new_L2P_k_modified);
			update_other_M2Ls_at_b(j, k, b, new_P2M_transpose_b_modified);

			if(j>2) {
				update_parent_L2P(j, k, new_L2P_k_modified);
				update_parent_P2M(j, b, new_P2M_transpose_b_modified);
			}

			tree[j][k].L2P = new_L2P_k_modified;
			tree[j][b].P2M = new_P2M_transpose_b_modified.adjoint();

			tree[j][k].NumLocals = tree[j][k].L2P.cols();
			tree[j][b].NumMultipoles = tree[j][b].P2M.rows();
			//////////////////////////////////////////////////////////////////////
			tree[j][b].M2L[k] = M2L_b_k.block(0,0,new_L2P_b_modified.cols(),new_P2M_transpose_k_modified.cols());
			// if (j==3 && b==6 && k==10) {
			// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree[j][k].M2L[b]: " << std::endl << tree[j][b].M2L[k] << std::endl << std::endl;
			// }
			update_other_M2Ls_at_k(j, b, k, new_L2P_b_modified);
			update_other_M2Ls_at_b(j, b, k, new_P2M_transpose_k_modified);

			if(j>2) {
				update_parent_L2P(j, b, new_L2P_b_modified);
				update_parent_P2M(j, k, new_P2M_transpose_k_modified);
			}

			tree[j][b].L2P = new_L2P_b_modified;
			tree[j][k].P2M = new_P2M_transpose_k_modified.adjoint();
			// std::cout << "j: " << j << "	k: " << k << "	P2P;	L2P.N: " << tree[j][k].L2P.norm() << "	P2M.N: " << tree[j][k].P2M.norm() << std::endl;
			// std::cout << "UtU: " << std::endl << tree[j][k].L2P.adjoint()*tree[j][k].L2P << std::endl << std::endl;
			// std::cout << "VVt: " << std::endl << tree[j][k].P2M*tree[j][k].P2M.adjoint() << std::endl << std::endl;
			// std::cout << "j: " << j << "	b: " << b << "	P2P;	L2P.N: " << tree[j][b].L2P.norm() << "	P2M.N: " << tree[j][b].P2M.norm() << std::endl;
			// if (j==3 && b==23) {
			// 	std::cout << "Vt:" << std::endl << new_P2M_transpose_b_modified << std::endl;
			// 	std::cout << "UtU: " << std::endl << tree[j][b].L2P.adjoint()*tree[j][b].L2P << std::endl << std::endl;
			// 	std::cout << "VVt: " << std::endl << tree[j][b].P2M*tree[j][b].P2M.adjoint() << std::endl << std::endl;
			// 	std::cout << "VtV: " << std::endl << tree[j][b].P2M.adjoint()*tree[j][b].P2M << std::endl << std::endl;
			// }

			tree[j][b].NumLocals = tree[j][b].L2P.cols();
			tree[j][k].NumMultipoles = tree[j][k].P2M.rows();
		}

		void compress_M2P_P2L(int j, int k, int b) {
			Mat new_L2P_k;
			Mat new_L2P_k_modified;
			Mat M2L_k_b;
			// std::cout << "here" << std::endl;
			compress_M2P(j, k, b, new_L2P_k, M2L_k_b);
			// std::cout << "compress_M2P done" << std::endl;
			Mat new_P2M_transpose_k;
			Mat new_P2M_transpose_k_modified;
			Mat M2L_b_k;
			compress_P2L(j, b, k, new_P2M_transpose_k, M2L_b_k);
			// std::cout << "compress_P2L done" << std::endl;

			if (new_L2P_k.cols() < new_P2M_transpose_k.cols()) {
				new_P2M_transpose_k_modified = new_P2M_transpose_k.block(0,0,new_P2M_transpose_k.rows(),new_L2P_k.cols());
				new_L2P_k_modified = new_L2P_k;
			}
			else {
				new_L2P_k_modified = new_L2P_k.block(0,0,new_L2P_k.rows(),new_P2M_transpose_k.cols());
				new_P2M_transpose_k_modified = new_P2M_transpose_k;
			}

			// std::cout << "L2P and P2M cutting done" << std::endl;

			tree[j][k].M2L[b] = M2L_k_b.block(0,0,new_L2P_k_modified.cols(),M2L_k_b.cols());
			// if (j==3 && k==6 && b==10) {
			// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree[j][k].M2L[b]: " << std::endl << tree[j][k].M2L[b] << std::endl << std::endl;
			// }
			update_other_M2Ls_at_k(j, k, b, new_L2P_k_modified);
			// std::cout << "update_other_M2Ls_at_k done" << std::endl;
			if(j>2) {
				update_parent_L2P(j, k, new_L2P_k_modified);
			}
			// std::cout << "update_parent_L2P done" << std::endl;
			tree[j][k].L2P = new_L2P_k_modified;
			tree[j][k].NumLocals = tree[j][k].L2P.cols();

			///////////////////////////////////////////////////////
			tree[j][b].M2L[k] = M2L_b_k.block(0,0,M2L_b_k.rows(),new_P2M_transpose_k_modified.cols());
			// if (j==3 && b==6 && k==10) {
			// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree[j][k].M2L[b]: " << std::endl << tree[j][b].M2L[k] << std::endl << std::endl;
			// }
			update_other_M2Ls_at_b(j, b, k, new_P2M_transpose_k_modified);
			// std::cout << "update_other_M2Ls_at_b done" << std::endl;
			if(j>2) {
				update_parent_P2M(j, k, new_P2M_transpose_k_modified);
			}
			// std::cout << "update_parent_P2M done" << std::endl;
			tree[j][k].P2M = new_P2M_transpose_k_modified.adjoint();
			// std::cout << "j: " << j << "	k: " << k << "	M2P;	L2P.N: " << tree[j][k].L2P.norm() << "	P2M.N: " << tree[j][k].P2M.norm() << std::endl;
			// std::cout << "UtU: " << std::endl << tree[j][k].L2P.adjoint()*tree[j][k].L2P << std::endl << std::endl;
			// std::cout << "VVt: " << std::endl << tree[j][k].P2M*tree[j][k].P2M.adjoint() << std::endl << std::endl;

			tree[j][k].NumMultipoles = tree[j][k].P2M.rows();
		}

		void compress_P2P(int j, int k, int b, Mat &new_L2P, Mat &new_P2M_transpose, Mat &M2L) {
			Mat A_col(tree[j][k].L2P.rows(), tree[j][k].L2P.cols()+tree[j][k].P2P[b].cols());
			A_col << tree[j][k].L2P, tree[j][k].P2P[b];
			Eigen::ColPivHouseholderQR<Mat> qr_A_col(A_col.rows(), A_col.cols());
			qr_A_col.setThreshold(RRQR_threshold);
			qr_A_col.compute(A_col);
			Mat new_L2P_old = qr_A_col.householderQ() ; //new tree[j][k].L2P
			Mat R_col = qr_A_col.matrixR().triangularView<Eigen::Upper>();
			new_L2P = new_L2P_old.block(0,0,new_L2P_old.rows(),qr_A_col.rank());
			// R_col = R_col.block(0,0,qr_A_col.rank(),R_col.cols());
			R_col = new_L2P.adjoint()*A_col;

			Mat Rk_prime = new_L2P.adjoint()*tree[j][k].P2P[b];
			Mat Rk = new_L2P.adjoint()*tree[j][k].L2P;

			Mat A_row(tree[j][b].P2M.cols(), tree[j][b].P2M.rows()+Rk_prime.rows());
			A_row << tree[j][b].P2M.adjoint(), Rk_prime.adjoint();
			Eigen::ColPivHouseholderQR<Mat> qr_A_row(A_row.rows(), A_row.cols());
			qr_A_row.setThreshold(RRQR_threshold);
			qr_A_row.compute(A_row);
			Mat new_P2M_transpose_old = qr_A_row.householderQ(); //new tree[j][n].L2P
			Mat R_row = qr_A_row.matrixQR().triangularView<Eigen::Upper>();
			new_P2M_transpose = new_P2M_transpose_old.block(0,0,new_P2M_transpose_old.rows(),qr_A_row.rank());
			// R_row = R_row.block(0,0,qr_A_row.rank(),R_row.cols());
			R_row = new_P2M_transpose.adjoint()*A_row;

			// if (j==3 && k==19 && b==23) {
			// 	std::cout << "Vt:" << std::endl << new_P2M_transpose << std::endl;
			// }

			Mat Rn = new_P2M_transpose.adjoint()*tree[j][b].P2M.adjoint();
			Mat Rn_prime = new_P2M_transpose.adjoint()*Rk_prime.adjoint();
			M2L = Rk*tree[j][k].M2L[b]*Rn.adjoint() + Rn_prime.adjoint();
		}

		void compress_M2P(int j, int k, int b, Mat &new_L2P, Mat &M2L) {
			Mat A_col(tree[j][k].L2P.rows(), tree[j][k].L2P.cols()+tree[j][k].M2P[b].cols());
			// std::cout << "M2P here" << std::endl;
			// std::cout << "L2P.S: " << tree[j][k].L2P.rows() << ", " << tree[j][k].L2P.cols() << std::endl;
			// std::cout << "M2P.S: " << tree[j][k].M2P[b].rows() << ", " << tree[j][k].M2P[b].cols() << std::endl;

			A_col << tree[j][k].L2P, tree[j][k].M2P[b];
			// std::cout << "M2P here" << std::endl;
			Eigen::ColPivHouseholderQR<Mat> qr_A_col(A_col.rows(), A_col.cols());
			qr_A_col.setThreshold(RRQR_threshold);
			qr_A_col.compute(A_col);
			Mat new_L2P_old = qr_A_col.householderQ() ; //new tree[j][k].L2P
			Mat R_col = qr_A_col.matrixR().triangularView<Eigen::Upper>();
			// std::cout << "M2P here" << std::endl;
			new_L2P = new_L2P_old.block(0,0,new_L2P_old.rows(),qr_A_col.rank());
			// R_col = R_col.block(0,0,qr_A_col.rank(),R_col.cols());
			R_col = new_L2P.adjoint()*A_col;
			// std::cout << "M2P here" << std::endl;

			Mat Rk_prime = new_L2P.adjoint()*tree[j][k].M2P[b];
			Mat Rk = new_L2P.adjoint()*tree[j][k].L2P;
			M2L = Rk*tree[j][k].M2L[b] + Rk_prime;
		}

		void compress_P2L(int j, int k, int b, Mat &new_P2M_transpose, Mat &M2L) {
			Mat A_row(tree[j][b].P2M.cols(), tree[j][b].P2M.rows()+tree[j][k].P2L[b].rows());
			A_row << tree[j][b].P2M.adjoint(), tree[j][k].P2L[b].adjoint();
			Eigen::ColPivHouseholderQR<Mat> qr_A_row(A_row.rows(), A_row.cols());
			qr_A_row.setThreshold(RRQR_threshold);
			qr_A_row.compute(A_row);
			Mat new_P2M_transpose_old = qr_A_row.householderQ(); //new tree[j][b].L2P
			Mat R_row = qr_A_row.matrixQR().triangularView<Eigen::Upper>();
			new_P2M_transpose = new_P2M_transpose_old.block(0,0,new_P2M_transpose_old.rows(),qr_A_row.rank());
			// R_row = R_row.block(0,0,qr_A_row.rank(),R_row.cols());
			R_row = new_P2M_transpose.adjoint()*A_row;

			Mat Rn = new_P2M_transpose.adjoint()*tree[j][b].P2M.adjoint();
			Mat Rn_prime = new_P2M_transpose.adjoint()*tree[j][k].P2L[b].adjoint();
			M2L = tree[j][k].M2L[b]*Rn.adjoint() + Rn_prime.adjoint();
		}


		void update_other_M2Ls_at_k(int j, int k, int b, Mat &new_L2P) {
			#pragma omp parallel for
			for(int in=0; in<16; ++in) {
				if(tree[j][k].innerNumbers[in] != -1) {
					int kIL = tree[j][k].innerNumbers[in];
					if(kIL == b) { //worked out in update_M2L_and_basis method
						continue;
					}
					tree[j][k].M2L[kIL] = new_L2P.adjoint() * tree[j][k].L2P * tree[j][k].M2L[kIL];
					// tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose(); // compressing P2P_{nn,k}
				}
			}
			#pragma omp parallel for
			for(int on=0; on<24; ++on) {
				if(tree[j][k].outerNumbers[on] != -1) {
					int kIL = tree[j][k].outerNumbers[on];
					if(kIL == b) { //worked out in update_M2L_and_basis method
						continue;
					}
					tree[j][k].M2L[kIL] = new_L2P.adjoint() * tree[j][k].L2P * tree[j][k].M2L[kIL];
					// tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose(); // compressing P2P_{nn,k}
				}
			}
			//k's Neighbors
			#pragma omp parallel for
			for(int n=0; n<8; ++n) {
				int kN = tree[j][k].neighborNumbers[n];
				if(kN != -1 && tree[j][kN].Eliminated && tree[j][k].Eliminated) {
					tree[j][k].M2L[kN] = new_L2P.adjoint() * tree[j][k].L2P * tree[j][k].M2L[kN];
					// tree[j][kN].M2L[k] = tree[j][k].M2L[kN].transpose(); // compressing P2P_{nn,k}
				}
			}
			//k's self
			if(tree[j][k].Eliminated) {
				tree[j][k].M2L[k] = new_L2P.adjoint() * tree[j][k].L2P * tree[j][k].M2L[k];
			}
		}

		void update_other_M2Ls_at_b(int j, int k, int b, Mat &new_P2M_transpose) {
			#pragma omp parallel for
			for(int in=0; in<16; ++in) {
				if(tree[j][b].innerNumbers[in] != -1) {
					int kIL = tree[j][b].innerNumbers[in];
					if(kIL == k) { //worked out in update_M2L_and_basis method
						continue;
					}
					tree[j][kIL].M2L[b] = tree[j][kIL].M2L[b] * tree[j][b].P2M * new_P2M_transpose;
					// if (j==3 && kIL==6 && b==10) {
					// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree[j][k].M2L[b]: " << std::endl << tree[j][kIL].M2L[b] << std::endl << std::endl;
					// }
					// tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose(); // compressing P2P_{nn,k}
				}
			}
			// std::cout << "inner done" << std::endl;
			#pragma omp parallel for
			for(int on=0; on<24; ++on) {
				if(tree[j][b].outerNumbers[on] != -1) {
					int kIL = tree[j][b].outerNumbers[on];
					if(kIL == b) { //worked out in update_M2L_and_basis method
						continue;
					}
					tree[j][kIL].M2L[b] = tree[j][kIL].M2L[b] * tree[j][b].P2M * new_P2M_transpose;
					// if (j==3 && kIL==6 && b==10) {
					// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree[j][k].M2L[b]: " << std::endl << tree[j][kIL].M2L[b] << std::endl << std::endl;
					// }
					// tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose(); // compressing P2P_{nn,k}
				}
			}
			// std::cout << "outer done" << std::endl;
			//k's Neighbors
			#pragma omp parallel for
			for(int n=0; n<8; ++n) {
				int kN = tree[j][b].neighborNumbers[n];
				if(kN != -1 && tree[j][kN].Eliminated && tree[j][b].Eliminated) {
					tree[j][kN].M2L[b] = tree[j][kN].M2L[b] * tree[j][b].P2M * new_P2M_transpose;
					// tree[j][kN].M2L[k] = tree[j][k].M2L[kN].transpose(); // compressing P2P_{nn,k}
				}
			}
			// std::cout << "neighbor done" << std::endl;
			//k's self
			if(tree[j][b].Eliminated) {
				tree[j][b].M2L[b] = tree[j][b].M2L[b] * tree[j][b].P2M * new_P2M_transpose;
			}
			// std::cout << "self done" << std::endl;
		}

		void update_parent_L2P(int j, int k, Mat &new_L2P) {
			int kP = k/4;//k's parent
			int c_k = k%4; // k is c_k^{th} child of k_parent

			Mat U[4];
			int row_index = 0;
			for (size_t c = 0; c < 4; c++) {
				U[c] = tree[j-1][kP].L2P.block(row_index,0,tree[j][kP*4+c].NumLocals,tree[j-1][kP].L2P.cols());
				row_index += tree[j][kP*4+c].NumLocals;
			}
			U[c_k] = new_L2P.adjoint() * tree[j][k].L2P * U[c_k];
			tree[j-1][kP].L2P = Mat(U[0].rows()+U[1].rows()+U[2].rows()+U[3].rows(), U[0].cols());
			tree[j-1][kP].L2P << U[0], U[1], U[2], U[3];
		}

		void update_parent_P2M(int j, int b, Mat &new_P2M_transpose) {
			int bP = b/4;//k's parent
			int c_b = b%4; // k is c_b^{th} child of k_parent

			Mat V[4];
			int col_index = 0;
			for (size_t c = 0; c < 4; c++) {
				V[c] = tree[j-1][bP].P2M.block(0,col_index,tree[j-1][bP].P2M.rows(),tree[j][bP*4+c].NumMultipoles);
				col_index += tree[j][bP*4+c].NumMultipoles;
			}
			V[c_b] = V[c_b] * tree[j][b].P2M * new_P2M_transpose;
			tree[j-1][bP].P2M = Mat(V[0].rows(), V[0].cols()+V[1].cols()+V[2].cols()+V[3].cols());
			tree[j-1][bP].P2M << V[0], V[1], V[2], V[3];
		}



	void initialise_rhs() {
		// #pragma omp parallel for
		for (size_t j = 2; j <= nLevels; j++) {
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				tree[j][k].local_rhs = Vec::Zero(tree[j][k].NumMultipoles);
				tree[j][k].multipole_rhs = Vec::Zero(tree[j][k].NumLocals);
			}
		}
	}

	void initialise_phase() {
		// #pragma omp parallel for
		for (size_t j = 2; j <= nLevels; j++) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				tree[j][k].NumMultipoles = tree[j][k].P2M.rows();
				tree[j][k].NumLocals = tree[j][k].L2P.cols();
        if (tree[j][k].NumMultipoles != tree[j][k].NumLocals) {
					// std::cout << "j: " << j << "	k: " << k << "	M: " << tree[j][k].NumMultipoles << "	L: " << tree[j][k].NumLocals << std::endl;
				}
				// std::cout << "j: " << j << "	k: " << k << "	M: " << tree[j][k].NumMultipoles << "	L: " << tree[j][k].NumLocals << std::endl;
			}
		}
		initialise_rhs();//initialise local_rhs and multipole_rhs
	}

	void initialise_P2P_Leaf_Level() {
		#pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
			// #pragma omp parallel for
			for(int n=0; n<8; ++n) {
				int nn = tree[nLevels][k].neighborNumbers[n];//P2P_{k,nn}
				if(nn != -1) {
					tree[nLevels][k].P2P[nn] = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][nn].chargeLocations);
				}
			}
			tree[nLevels][k].P2P[k] = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][k].chargeLocations);
		}
	}

	void initialise_P2P_NonLeafLevel(int j) {
		#pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			// #pragma omp parallel for
			for(int n=0; n<8; ++n) {
				int nn = tree[j][k].neighborNumbers[n];//P2P_{k,nn}
				if(nn != -1) {
					int n_rows = 0;
					int n_cols = 0;
					for(int a=0; a<4; ++a) {
						n_rows += tree[j+1][4*k+a].NumLocals;//tree[j+1][4*k+a].M2L[4*nn].rows();
					}
					for(int b=0; b<4; ++b) {
						n_cols += tree[j+1][4*nn+b].NumMultipoles;//tree[j+1][4*k].M2L[4*nn+b].cols();
					}
					tree[j][k].P2P[nn] = Mat::Zero(n_rows, n_cols);
					// if (j==3 && k==27 && nn==37) {
					// 	std::cout << "tree[j][k].P2P[nn].S: " << tree[j][k].P2P[nn].rows() << ", " << tree[j][k].P2P[nn].cols() << std::endl;
					// }
					int row_index = 0;
					int col_index = 0;
					for(int a=0; a<4; ++a) {
						for(int b=0; b<4; ++b) {
							tree[j][k].P2P[nn].block(row_index,col_index,tree[j+1][4*k+a].NumLocals,tree[j+1][4*nn+b].NumMultipoles) = tree[j+1][4*k+a].M2L[4*nn+b];
							// if(k==1 && nn==2) {
							// 	std::cout << "M2L: " << j+1 << ", " << 4*k+a << ", " << 4*nn+b << "----------------" << std::endl;
							// 	std::cout << tree[j+1][4*k+a].M2L[4*nn+b] << std::endl << std::endl;
							// }
							col_index += tree[j+1][4*nn+b].NumMultipoles;//tree[j+1][4*k+a].M2L[4*nn+b].cols();
						}
						row_index += tree[j+1][4*k+a].NumLocals;//tree[j+1][4*k+a].M2L[4*nn].rows();
						col_index = 0;
					}
					// std::cout << "k: " << k << "	nn: " << nn << "	P2P.n: " << tree[j][k].P2P[nn].norm() << std::endl;
					// if(k==1 && nn==2) {
					// 	std::cout << "P2P: -------------------------" << std::endl << tree[j][k].P2P[nn] << std::endl;
					// }
				}
			}
			// P2P_self
			{
				int nn = k;
				int n_rows = 0;
				int n_cols = 0;
				for(int a=0; a<4; ++a) {
					n_rows += tree[j+1][4*k+a].NumLocals;//tree[j+1][4*k+a].M2L[4*nn].rows();
				}
				for(int b=0; b<4; ++b) {
					n_cols += tree[j+1][4*nn+b].NumMultipoles;//tree[j+1][4*k].M2L[4*nn+b].cols();
				}
				tree[j][k].P2P[nn] = Mat::Zero(n_rows, n_cols);
				int row_index = 0;
				int col_index = 0;
				for(int a=0; a<4; ++a) {
					for(int b=0; b<4; ++b) {
						tree[j][k].P2P[nn].block(row_index,col_index,tree[j+1][4*k+a].NumLocals,tree[j+1][4*nn+b].NumMultipoles) = tree[j+1][4*k+a].M2L[4*nn+b];
						col_index += tree[j+1][4*nn+b].NumMultipoles;//tree[j+1][4*k+a].M2L[4*nn+b].cols();
					}
					row_index += tree[j+1][4*k+a].NumLocals;//tree[j+1][4*k+a].M2L[4*nn].rows();
					col_index = 0;
				}
			}
		}
    // std::cout << "-----------DONE-------------------------" << std::endl;
		// assign_particle_rhs_NonLeafLevel(j);
    // std::cout << "-----------DONE-------------------------" << std::endl;
	}

	void assign_particle_rhs_NonLeafLevel(int j) {
    #pragma omp parallel for
    for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
			int particle_rhs_length = 0;
			for (size_t c = 0; c < 4; c++) {
				// particle_rhs_length += tree[j+1][4*k+c].NumLocals;
        particle_rhs_length += tree[j+1][4*k+c].multipole_rhs.size();
        if(k==10) {
          // std::cout << "c: " << c <<  " :" << tree[j+1][4*k+0].multipole_rhs.size() << std::endl;
        }
			}
			tree[j][k].particle_rhs = Vec(particle_rhs_length);
			tree[j][k].particle_rhs << tree[j+1][4*k+0].multipole_rhs, tree[j+1][4*k+1].multipole_rhs, tree[j+1][4*k+2].multipole_rhs, tree[j+1][4*k+3].multipole_rhs;
		}
	}

  void rhs_eliminate_cluster(int j, int k) {
		rhs_eliminate_x(j, k);
		// std::cout << "eliminate_x done" << std::endl;
		// if (tree[j][k].P2M.rows() > 0) {
			rhs_eliminate_z(j, k);
			// std::cout << "eliminate_z done" << std::endl;
		// }
		// else {
		// 	std::cout << "eliminate_z NOT done" << std::endl;
		// }
	}

	void eliminate_cluster(int j, int k) {
		eliminate_x(j, k);
		// std::cout << "eliminate_x done" << std::endl;
		// if (tree[j][k].P2M.rows() > 0) {
			eliminate_z(j, k);
			// std::cout << "eliminate_z done" << std::endl;
		// }
		// else {
		// 	std::cout << "eliminate_z NOT done" << std::endl;
		// }
	}

	void eliminate_x(int j, int k) {
		// in P2P(j,k) we have terms: P2P_SELF, U, P2P_neighbor;
		// self is eliminated, so we have fill-ins/updates with respect to U, P2P_neighbor
		// in equation U.transpose()(j,k)
		// if (j==3 && k==39) {
			// std::cout << "tree[j][k].P2P[k].S: " << tree[j][k].P2P[k].rows() << ", " << tree[j][k].P2P[k].cols() << std::endl;
		// }
		if (tree[j][k].P2P[k].size() != 0) {
			Eigen::PartialPivLU<Mat> P2P_self_QR = tree[j][k].P2P[k].lu();
			filings_in_equation_P2M_due_to_x(j,k,P2P_self_QR);//equation z
			// std::cout << "filings_in_equation_P2M_due_to_x done" << std::endl;
			// std::cout << "here" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < 8; n++) { //loop over P2P(j,N(k))
				int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
				if(nn != -1) {
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_x START" << std::endl;
					filings_in_equation_P2P_due_to_x(j,k,nn,P2P_self_QR);
				}
			}
		}
		else {
			// std::cout << "here" << std::endl;
			// Eigen::ColPivHouseholderQR<Mat> P2P_self_QR = tree[j][k].P2P[k].colPivHouseholderQr();
			// std::cout << "here" << std::endl;
			filings_in_equation_P2M_due_to_x(j,k);//equation z
			// std::cout << "filings_in_equation_P2M_due_to_x done" << std::endl;
			// std::cout << "here" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < 8; n++) { //loop over P2P(j,N(k))
				int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
				if(nn != -1) {
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_x START" << std::endl;
					filings_in_equation_P2P_due_to_x(j,k,nn);
				}
			}
		}
	}

  void rhs_eliminate_x(int j, int k) {
		// in P2P(j,k) we have terms: P2P_SELF, U, P2P_neighbor;
		// self is eliminated, so we have fill-ins/updates with respect to U, P2P_neighbor
		// in equation U.transpose()(j,k)
		// if (j==3 && k==39) {
			// std::cout << "tree[j][k].P2P[k].S: " << tree[j][k].P2P[k].rows() << ", " << tree[j][k].P2P[k].cols() << std::endl;
		// }
		if (tree[j][k].P2P[k].size() != 0) {
      // std::cout << "x here !=0" << std::endl;
			Eigen::PartialPivLU<Mat> P2P_self_QR = tree[j][k].P2P[k].lu();
			rhs_filings_in_equation_P2M_due_to_x(j,k,P2P_self_QR);//equation z
      // std::cout << "x here" << std::endl;
      // std::cout << "filings_in_equation_P2M_due_to_x done" << std::endl;
			// std::cout << "here" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < 8; n++) { //loop over P2P(j,N(k))
				int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
				if(nn != -1) {
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_x START" << std::endl;
					rhs_filings_in_equation_P2P_due_to_x(j,k,nn,P2P_self_QR);
          // std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_x DONE" << std::endl;
				}
			}
		}
		else {
			// std::cout << "here 0" << std::endl;
			// Eigen::ColPivHouseholderQR<Mat> P2P_self_QR = tree[j][k].P2P[k].colPivHouseholderQr();
			// std::cout << "here" << std::endl;
			rhs_filings_in_equation_P2M_due_to_x(j,k);//equation z
			// std::cout << "filings_in_equation_P2M_due_to_x done" << std::endl;
			// std::cout << "here" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < 8; n++) { //loop over P2P(j,N(k))
				int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
				if(nn != -1) {
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_x START" << std::endl;
					rhs_filings_in_equation_P2P_due_to_x(j,k,nn);
				}
			}
		}
	}

////////////////////////////////////////////////////////////////////
void filings_in_equation_P2M_due_to_x(int j, int k) {
	filing_in_equation_P2M_due_to_L2P(j,k);//U fill-in/update; in P2P(j,n) equation
	// std::cout << "filing_in_equation_P2M_due_to_L2P done" << std::endl;
	filing_in_equation_P2M_due_to_P2P_neighbors(j,k);//P2P_neighbor fill-in/update; in P2P(j,n) equation
	// std::cout << "filing_in_equation_P2M_due_to_P2P_neighbors done" << std::endl;
	// rhs_update_in_equation_P2M_due_to_x(j,k);
	// std::cout << "rhs_update_in_equation_P2M_due_to_x done" << std::endl;
}

void rhs_filings_in_equation_P2M_due_to_x(int j, int k) {
	rhs_update_in_equation_P2M_due_to_x(j,k);
	// std::cout << "rhs_update_in_equation_P2M_due_to_x done" << std::endl;
}

void rhs_update_in_equation_P2M_due_to_x(int j, int k) {
	tree[j][k].local_rhs = Vec(tree[j][k].P2M.rows());//- tree[j][k].P2M * P2P_self_QR.solve(tree[j][k].particle_rhs);
}

void filing_in_equation_P2M_due_to_L2P(int j, int k) {
	tree[j][k].L2M_f = Mat(tree[j][k].P2M.rows(),tree[j][k].L2P.cols());//-tree[j][k].P2M * P2P_self_QR.solve(tree[j][k].L2P);
}

void filing_in_equation_P2M_due_to_P2P_neighbors(int j, int k) {
	for (size_t n = 0; n < 8; n++) {//neighbors of k
		int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
		if (nn != -1) {
			if (!tree[j][nn].Eliminated) {
				tree[j][k].P2M_f[nn] = Mat(tree[j][k].P2M.rows(),tree[j][k].P2P[nn].cols());//-tree[j][k].P2M * P2P_self_QR.solve(tree[j][k].P2P[nn]);
			}
			else {
				tree[j][k].M2M_f[nn] = Mat(tree[j][k].P2M.rows(),tree[j][k].M2P[nn].cols());//-tree[j][k].P2M * P2P_self_QR.solve(tree[j][k].M2P[nn]);
			}
		}
	}
}

void filings_in_equation_P2P_due_to_x(int j, int k, int nn) {
	filing_due_to_L2P(j,k,nn);//U fill-in/update; in P2P(j,n) equation
	// std::cout << "filing_due_to_L2P done" << std::endl;
	filing_due_to_P2P_neighbors(j,k,nn);//P2P_neighbor fill-in/update; in P2P(j,n) equation
	// std::cout << "filing_due_to_P2P_neighbors done" << std::endl;
	// rhs_update_in_equation_P2P_due_to_x(j,k,nn);
	// std::cout << "rhs_update_in_equation_P2P_due_to_x done" << std::endl;
}

void rhs_filings_in_equation_P2P_due_to_x(int j, int k, int nn) {
	rhs_update_in_equation_P2P_due_to_x(j,k,nn);
	// std::cout << "rhs_update_in_equation_P2P_due_to_x done" << std::endl;
}

void rhs_update_in_equation_P2P_due_to_x(int j, int k, int nn) {
	if (!tree[j][nn].Eliminated) {
		tree[j][nn].particle_rhs = tree[j][nn].particle_rhs;// - tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].particle_rhs);
	}
	else {
		tree[j][nn].multipole_rhs = tree[j][nn].multipole_rhs;// - tree[j][nn].P2L[k] * P2P_self_QR.solve(tree[j][k].particle_rhs);
	}
}

void filing_due_to_L2P(int j, int k, int nn) {//in P2P(j,nn) equation; nn:eq number
	if (!tree[j][nn].Eliminated) {
		tree[j][nn].L2P_f[k] = Mat(tree[j][nn].P2P[k].rows(),tree[j][k].L2P.cols());//-tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].L2P); //temporary; this gets eliminated when z is eliminated
	}
	else {
		// if (k==1 && nn==0) {
		// 	std::cout << "tree[j][nn].P2L[k].S: " << tree[j][nn].P2L[k].rows() << ", " << tree[j][nn].P2L[k].cols() << std::endl;
		// 	std::cout << "tree[j][k].P2P[k].S: " << tree[j][k].P2P[k].rows() << ", " << tree[j][k].P2P[k].cols() << std::endl;
		// 	std::cout << "tree[j][k].L2P.S: " << tree[j][nn].L2P.rows() << ", " << tree[j][nn].L2P.cols() << std::endl;
		// }
		tree[j][nn].L2P_f[k] = Mat(tree[j][nn].P2L[k].rows(),tree[j][k].L2P.cols());//-tree[j][nn].P2L[k] * P2P_self_QR.solve(tree[j][k].L2P);
	}// 0-K21 K11^{-1} U1
}

void filing_due_to_P2P_neighbors(int j, int k, int nn) {//in P2P(j,nn) equation
	// Eq: P2P_n
	// workout for fill-ins due to all neighbors of (j,k)
	for (size_t p = 0; p < 8; p++) { //loop over all neighbors of (j,k)
		int pp = tree[j][k].neighborNumbers[p];
		if (!tree[j][nn].Eliminated) {
			if (pp != -1) {
				if (!tree[j][pp].Eliminated) {
					if (tree[j][nn].P2P[pp].size() == 0) {
						tree[j][nn].P2P[pp] = Mat::Zero(tree[j][nn].P2P[k].rows(),tree[j][k].P2P[pp].cols());//-tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].P2P[pp]);
					}
					else {
						tree[j][nn].P2P[pp] = tree[j][nn].P2P[pp];// - tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].P2P[pp]);
					}
					if(is_well_separated(j,nn,pp)) {
						if (nn < pp) { // using symmetry
							if(tree[j][nn].P2P[pp].norm() > 1e-10) {
								std::pair<int, int> g(nn,pp);
								P2P.push_back(g);
							}
						}
					}
				}
				else {
					if (tree[j][nn].M2P[pp].size() == 0) {
						tree[j][nn].M2P[pp] = Mat(tree[j][nn].P2P[k].rows(),tree[j][k].M2P[pp].cols());//-tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].M2P[pp]);
					}
					else {
						tree[j][nn].M2P[pp] = tree[j][nn].M2P[pp];// - tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].M2P[pp]);
					}
				}
			}
		}
		else {
			if (pp != -1) {
				if (!tree[j][pp].Eliminated) {
					if (tree[j][nn].P2L[pp].size() == 0) {
						tree[j][nn].P2L[pp] = Mat(tree[j][nn].P2L[k].rows(),tree[j][k].P2P[pp].cols());//-tree[j][nn].P2L[k] * P2P_self_QR.solve(tree[j][k].P2P[pp]);
					}
					else {
						tree[j][nn].P2L[pp] = tree[j][nn].P2L[pp];// - tree[j][nn].P2L[k] * P2P_self_QR.solve(tree[j][k].P2P[pp]);
					}
					if(is_well_separated(j,nn,pp)) {
						if (tree[j][nn].P2L[pp].norm() > 1e-10) {
							std::pair<int, int> g(nn,pp);
							P2L_M2P.push_back(g);
						}
					}
				}
				else {
					tree[j][nn].M2L[pp] = tree[j][nn].M2L[pp];// - tree[j][nn].P2L[k] * P2P_self_QR.solve(tree[j][k].M2P[pp]);
				}
			}
		}
		// p : 2,4; nn: 2
		// Eq 2:(n=2)
		// K22-K21 K11^{-1}K12 :p=2; K22 update
		// -K21 K11^{-1}K14 :p=4; K24 update
		// Eq 4:(n=4)
		// -K41 K11^{-1}K12 :p=2		 -Kn1 K11^{-1}K1p
		// K44-K41 K11^{-1}K14 :p=4		Knp-Kn1 K11^{-1}K1p
	}
}


///////////////////////////////////////////////////////////////////
void rhs_filings_in_equation_P2M_due_to_x(int j, int k, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
  rhs_update_in_equation_P2M_due_to_x(j,k,P2P_self_QR);
  // std::cout << "rhs_update_in_equation_P2M_due_to_x done" << std::endl;
}

	void filings_in_equation_P2M_due_to_x(int j, int k, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		filing_in_equation_P2M_due_to_L2P(j,k,P2P_self_QR);//U fill-in/update; in P2P(j,n) equation
		// std::cout << "filing_in_equation_P2M_due_to_L2P done" << std::endl;
		filing_in_equation_P2M_due_to_P2P_neighbors(j,k,P2P_self_QR);//P2P_neighbor fill-in/update; in P2P(j,n) equation
		// std::cout << "filing_in_equation_P2M_due_to_P2P_neighbors done" << std::endl;
		// rhs_update_in_equation_P2M_due_to_x(j,k,P2P_self_QR);
		// std::cout << "rhs_update_in_equation_P2M_due_to_x done" << std::endl;
	}

	void rhs_update_in_equation_P2M_due_to_x(int j, int k, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		tree[j][k].local_rhs = - tree[j][k].P2M * P2P_self_QR.solve(tree[j][k].particle_rhs);
	}

	void filing_in_equation_P2M_due_to_L2P(int j, int k, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		tree[j][k].L2M_f = -tree[j][k].P2M * P2P_self_QR.solve(tree[j][k].L2P);
	}

	void filing_in_equation_P2M_due_to_P2P_neighbors(int j, int k, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		// if (j==3 && k==27) {
		// 	std::cout << "tttttt" << std::endl;
		// }
		for (size_t n = 0; n < 8; n++) {//neighbors of k
			int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
			// if (j==3 && k==27) {
			// 	std::cout << "k: " << k << "	nn: " << nn << std::endl;
			// 	std::cout << "is_well_separated: " << is_well_separated(j,k,nn) << std::endl;
			// }
			if (nn != -1) {
				if (!tree[j][nn].Eliminated) {
					// if (j==3 && k==27) {
					// 	std::cout << "tree[j][k].P2M.S: " << tree[j][k].P2M.rows() << ", " << tree[j][k].P2M.cols() << std::endl;
					// 	std::cout << "tree[j][k].P2P[k].S: " << tree[j][k].P2P[k].rows() << ", " << tree[j][k].P2P[k].cols() << std::endl;
					// 	std::cout << "tree[j][k].P2P[nn].S: " << tree[j][k].P2P[nn].rows() << ", " << tree[j][k].P2P[nn].cols() << std::endl;
					// }
					tree[j][k].P2M_f[nn] = -tree[j][k].P2M * P2P_self_QR.solve(tree[j][k].P2P[nn]);
				}
				else {
					// if (j==3 && k==27) {
					// 	std::cout << "tree[j][k].P2M.S: " << tree[j][k].P2M.rows() << ", " << tree[j][k].P2M.cols() << std::endl;
					// 	std::cout << "tree[j][k].P2P[k].S: " << tree[j][k].P2P[k].rows() << ", " << tree[j][k].P2P[k].cols() << std::endl;
					// 	std::cout << "tree[j][k].M2P[nn].S: " << tree[j][k].M2P[nn].rows() << ", " << tree[j][k].M2P[nn].cols() << std::endl;
					// }
					tree[j][k].M2M_f[nn] = -tree[j][k].P2M * P2P_self_QR.solve(tree[j][k].M2P[nn]);
				}
			}
		}
	}

	void filings_in_equation_P2P_due_to_x(int j, int k, int nn, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		filing_due_to_L2P(j,k,nn,P2P_self_QR);//U fill-in/update; in P2P(j,n) equation
		// std::cout << "filing_due_to_L2P done" << std::endl;
		filing_due_to_P2P_neighbors(j,k,nn,P2P_self_QR);//P2P_neighbor fill-in/update; in P2P(j,n) equation
		// std::cout << "filing_due_to_P2P_neighbors done" << std::endl;
		// rhs_update_in_equation_P2P_due_to_x(j,k,nn,P2P_self_QR);
		// std::cout << "rhs_update_in_equation_P2P_due_to_x done" << std::endl;
	}

  void rhs_filings_in_equation_P2P_due_to_x(int j, int k, int nn, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
		rhs_update_in_equation_P2P_due_to_x(j,k,nn,P2P_self_QR);
		// std::cout << "rhs_update_in_equation_P2P_due_to_x done" << std::endl;
	}

	void rhs_update_in_equation_P2P_due_to_x(int j, int k, int nn, Eigen::PartialPivLU<Mat>& P2P_self_QR) {
    // if (j==2 && k==8 && nn==10) {
    //   std::cout << "tree[j][nn].Eliminated: " << tree[j][nn].Eliminated << std::endl;
    //   std::cout << "tree[j][nn].P2P[k]: " << tree[j][nn].P2P[k].rows() << "," << tree[j][nn].P2P[k].cols() << std::endl;
    //   std::cout << "tree[j][nn].P2L[k]: " << tree[j][nn].P2L[k].rows() << "," << tree[j][nn].P2L[k].cols() << std::endl;
    //   std::cout << "tree[j][nn].particle_rhs: " << tree[j][nn].particle_rhs.size() << std::endl;
    //   std::cout << "tree[j][nn].multipole_rhs: " << tree[j][nn].multipole_rhs.size() << std::endl;
    // }
		if (!tree[j][nn].Eliminated) {
			tree[j][nn].particle_rhs = tree[j][nn].particle_rhs - tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].particle_rhs);
		}
		else {
			tree[j][nn].multipole_rhs = tree[j][nn].multipole_rhs - tree[j][nn].P2L[k] * P2P_self_QR.solve(tree[j][k].particle_rhs);
		}
	}

	void filing_due_to_L2P(int j, int k, int nn, Eigen::PartialPivLU<Mat>& P2P_self_QR) {//in P2P(j,nn) equation; nn:eq number
		if (!tree[j][nn].Eliminated) {
			tree[j][nn].L2P_f[k] = -tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].L2P); //temporary; this gets eliminated when z is eliminated
		}
		else {
			// if (k==1 && nn==0) {
			// 	std::cout << "tree[j][nn].P2L[k].S: " << tree[j][nn].P2L[k].rows() << ", " << tree[j][nn].P2L[k].cols() << std::endl;
			// 	std::cout << "tree[j][k].P2P[k].S: " << tree[j][k].P2P[k].rows() << ", " << tree[j][k].P2P[k].cols() << std::endl;
			// 	std::cout << "tree[j][k].L2P.S: " << tree[j][nn].L2P.rows() << ", " << tree[j][nn].L2P.cols() << std::endl;
			// }
			tree[j][nn].L2P_f[k] = -tree[j][nn].P2L[k] * P2P_self_QR.solve(tree[j][k].L2P);
		}// 0-K21 K11^{-1} U1
	}

	void filing_due_to_P2P_neighbors(int j, int k, int nn, Eigen::PartialPivLU<Mat>& P2P_self_QR) {//in P2P(j,nn) equation
		// Eq: P2P_n
		// workout for fill-ins due to all neighbors of (j,k)
		for (size_t p = 0; p < 8; p++) { //loop over all neighbors of (j,k)
			int pp = tree[j][k].neighborNumbers[p];
			// if (j==3 && k==39 && nn==36) {
			// 	std::cout << "pp: " << pp << "	E: " << tree[j][pp].Eliminated << std::endl;
			// 	std::cout << "nn: " << nn << "	E: " << tree[j][nn].Eliminated << std::endl;
			// 	std::cout << "is_well_separated(j,k,nn): " << is_well_separated(j,k,nn) << std::endl;
			// 	std::cout << "is_well_separated(j,nn,pp): " << is_well_separated(j,nn,pp) << std::endl;
			// 	std::cout << "is_well_separated(j,k,pp): " << is_well_separated(j,k,pp) << std::endl;
			// }
			if (!tree[j][nn].Eliminated) {
				if (pp != -1) {
					if (!tree[j][pp].Eliminated) {
						if (tree[j][nn].P2P[pp].size() == 0) {
							tree[j][nn].P2P[pp] = -tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].P2P[pp]);
						}
						else {
							tree[j][nn].P2P[pp] = tree[j][nn].P2P[pp] - tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].P2P[pp]);
						}
						if(is_well_separated(j,nn,pp)) {
							if (nn < pp) { // using symmetry
								if(tree[j][nn].P2P[pp].norm() > 1e-10) {
									std::pair<int, int> g(nn,pp);
									P2P.push_back(g);
								}
							}
						}
					}
					else {
						if (tree[j][nn].M2P[pp].size() == 0) {
							tree[j][nn].M2P[pp] = -tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].M2P[pp]);
						}
						else {
							tree[j][nn].M2P[pp] = tree[j][nn].M2P[pp] - tree[j][nn].P2P[k] * P2P_self_QR.solve(tree[j][k].M2P[pp]);
						}
					}
				}
			}
			else {
				if (pp != -1) {
					if (!tree[j][pp].Eliminated) {
						if (tree[j][nn].P2L[pp].size() == 0) {
							tree[j][nn].P2L[pp] = -tree[j][nn].P2L[k] * P2P_self_QR.solve(tree[j][k].P2P[pp]);
						}
						else {
							tree[j][nn].P2L[pp] = tree[j][nn].P2L[pp] - tree[j][nn].P2L[k] * P2P_self_QR.solve(tree[j][k].P2P[pp]);
						}
						if(is_well_separated(j,nn,pp)) {
							if (tree[j][nn].P2L[pp].norm() > 1e-10) {
								std::pair<int, int> g(nn,pp);
								P2L_M2P.push_back(g);
							}
						}
					}
					else {
						// if (j==3 && k==39 && nn==36 && pp==38) {
						// 	std::cout << "tree[j][nn].M2L[pp].S: " << tree[j][nn].M2L[pp].rows() << ", " << tree[j][nn].M2L[pp].cols() << std::endl;
						// 	std::cout << "tree[j][nn].P2L[k].S: " << tree[j][nn].P2L[k].rows() << ", " << tree[j][nn].P2L[k].cols() << std::endl;
						// 	std::cout << "tree[j][k].P2P[k].S: " << tree[j][k].P2P[k].rows() << ", " << tree[j][nn].P2P[k].cols() << std::endl;
						// 	std::cout << "tree[j][k].M2P[pp].S: " << tree[j][k].M2P[pp].rows() << ", " << tree[j][k].M2P[pp].cols() << std::endl;
						// }
						tree[j][nn].M2L[pp] = tree[j][nn].M2L[pp] - tree[j][nn].P2L[k] * P2P_self_QR.solve(tree[j][k].M2P[pp]);
						// if (j==3 && nn==6 && pp==10) {
						// 	std::cout << "eliminate_x: k: " << k << std::endl;
						// 	std::cout << "tree[j][nn].P2L[k]: " << std::endl << tree[j][nn].P2L[k] << std::endl << std::endl;
						// 	std::cout << "tree[j][k].P2P[k]: " << std::endl << tree[j][k].P2P[k] << std::endl << std::endl;
						// 	std::cout << "tree[j][k].M2P[pp]: " << std::endl << tree[j][k].M2P[pp] << std::endl << std::endl;
						// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree[j][nn].M2L[pp]: " << std::endl << tree[j][nn].M2L[pp] << std::endl << std::endl;
						// }
					}
				}
			}
			// p : 2,4; nn: 2
			// Eq 2:(n=2)
			// K22-K21 K11^{-1}K12 :p=2; K22 update
			// -K21 K11^{-1}K14 :p=4; K24 update
			// Eq 4:(n=4)
			// -K41 K11^{-1}K12 :p=2		 -Kn1 K11^{-1}K1p
			// K44-K41 K11^{-1}K14 :p=4		Knp-Kn1 K11^{-1}K1p
		}
	}

	bool is_well_separated(int j, int k, int l) {
		if(k==l) { //self
			return false;
		}
		for (size_t i = 0; i < 8; i++) {
			int nn = tree[j][k].neighborNumbers[i];
			if(l==nn) {
				return false;
			}
		}
		return true;
	}

  void eliminate_z(int j, int k) {
		// if(j==2 && k==0) {
		// 	std::cout << "tree[j][k].L2M_f.S: " << tree[j][k].L2M_f.rows() << ", " << tree[j][k].L2M_f.cols() << std::endl;
		// }
		if (tree[j][k].L2M_f.size() > 0) {
			Eigen::PartialPivLU<Mat> L2M_f_QR = tree[j][k].L2M_f.lu();
			// std::cout << "there" << std::endl;
			filings_in_equation_M2L_due_to_z(j,k,L2M_f_QR);
			// std::cout << "filings_in_equation_M2L_due_to_z done" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < 8; n++) { //loop over P2P(j,N(k))
				int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
				if(nn != -1) {
					filings_in_equation_P2P_due_to_z(j,k,nn,L2M_f_QR);
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_z done" << std::endl;
				}
			}
		}
		else {
			// Eigen::ColPivHouseholderQR<Mat> L2M_f_QR = tree[j][k].L2M_f.colPivHouseholderQr();
			// std::cout << "there" << std::endl;
			filings_in_equation_M2L_due_to_z(j,k);
			// std::cout << "filings_in_equation_M2L_due_to_z done" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < 8; n++) { //loop over P2P(j,N(k))
				int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
				if(nn != -1) {
					filings_in_equation_P2P_due_to_z(j,k,nn);
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_z done" << std::endl;
				}
			}
		}
	}

	void rhs_eliminate_z(int j, int k) {
		// if(j==2 && k==0) {
		// 	std::cout << "tree[j][k].L2M_f.S: " << tree[j][k].L2M_f.rows() << ", " << tree[j][k].L2M_f.cols() << std::endl;
		// }
		if (tree[j][k].L2M_f.size() > 0) {
			Eigen::PartialPivLU<Mat> L2M_f_QR = tree[j][k].L2M_f.lu();
			// std::cout << "THERE " << k << std::endl;
			rhs_filings_in_equation_M2L_due_to_z(j,k,L2M_f_QR);
			// std::cout << "filings_in_equation_M2L_due_to_z done" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < 8; n++) { //loop over P2P(j,N(k))
				int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
				if(nn != -1) {
          // std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_z START" << std::endl;
					rhs_filings_in_equation_P2P_due_to_z(j,k,nn,L2M_f_QR);
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_z done" << std::endl;
				}
			}
		}
		else {
			// Eigen::ColPivHouseholderQR<Mat> L2M_f_QR = tree[j][k].L2M_f.colPivHouseholderQr();
			// std::cout << "there" << std::endl;
			rhs_filings_in_equation_M2L_due_to_z(j,k);
			// std::cout << "filings_in_equation_M2L_due_to_z done" << std::endl;
			// in equation P2P(j,N(k))
			for (size_t n = 0; n < 8; n++) { //loop over P2P(j,N(k))
				int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
				if(nn != -1) {
					rhs_filings_in_equation_P2P_due_to_z(j,k,nn);
					// std::cout << "nn: " << nn << "	filings_in_equation_P2P_due_to_z done" << std::endl;
				}
			}
		}
	}

	///////////////////////////////////////////////
  void rhs_filings_in_equation_M2L_due_to_z(int j, int k) {
		rhs_update_in_equation_M2L(j,k);
		// std::cout << "rhs_update_in_equation_M2L done" << std::endl;
	}

	void filings_in_equation_M2L_due_to_z(int j, int k) {
		filings_in_equation_M2L_due_to_y(j,k);
		// std::cout << "filings_in_equation_M2L_due_to_y done" << std::endl;
		filings_in_equation_M2L_due_to_neighbors(j,k);
		// std::cout << "filings_in_equation_M2L_due_to_neighbors done" << std::endl;
		// rhs_update_in_equation_M2L(j,k);
		// std::cout << "rhs_update_in_equation_M2L done" << std::endl;
	}

	void rhs_update_in_equation_M2L(int j, int k) {
		tree[j][k].multipole_rhs = Vec(0);//L2M_f_QR.solve(tree[j][k].local_rhs);
	}

	void filings_in_equation_M2L_due_to_y(int j, int k) {
		tree[j][k].M2L[k] = Mat(0,tree[j][k].NumMultipoles);//L2M_f_QR.solve(-Mat::Identity(tree[j][k].NumMultipoles, tree[j][k].NumMultipoles));
	}

	void filings_in_equation_M2L_due_to_neighbors(int j, int k) {
		for (size_t n = 0; n < 8; n++) {
			int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
			if(nn != -1) {
				if (!tree[j][nn].Eliminated) {
					tree[j][k].P2L[nn] = Mat(0,tree[j][k].P2M_f[nn].cols());//L2M_f_QR.solve(tree[j][k].P2M_f[nn]);
					if(is_well_separated(j,k,nn)) {
						if (tree[j][k].P2L[nn].norm() > 1e-10) {
							// std::cout << "compress_P2L here" << std::endl;
							std::pair<int, int> g(k,nn);
							P2L_M2P.push_back(g);
							// compress_P2L(j,k,nn);
							// std::cout << "there" << std::endl;
							// std::cout << "j: " << j << "	k: " << k << "	nn: " << nn << "	compress_P2L" << std::endl;
						}
						// tree[j][k].P2L[nn] = Mat(0,0);
					}
				}
				else {
					tree[j][k].M2L[nn] = Mat(0,tree[j][k].M2M_f[nn].cols());//L2M_f_QR.solve(tree[j][k].M2M_f[nn]);
				}
			}
		}
	}

  void rhs_filings_in_equation_P2P_due_to_z(int j, int k, int nn) {
		rhs_update_in_equation_P2P_due_to_z(j,k,nn);
		// std::cout << "rhs_update_in_equation_P2P_due_to_z done" << std::endl;
	}

	void filings_in_equation_P2P_due_to_z(int j, int k, int nn) {
		filings_in_equation_P2P_due_to_y(j,k,nn);
		// std::cout << "filings_in_equation_P2P_due_to_y done" << std::endl;
		filings_in_equation_P2P_due_to_neighbors(j,k,nn);
		// std::cout << "filings_in_equation_P2P_due_to_neighbors done" << std::endl;
		// rhs_update_in_equation_P2P_due_to_z(j,k,nn);
		// std::cout << "rhs_update_in_equation_P2P_due_to_z done" << std::endl;
	}

	void rhs_update_in_equation_P2P_due_to_z(int j, int k, int nn) {
		if (!tree[j][nn].Eliminated) {
			tree[j][nn].particle_rhs = tree[j][nn].particle_rhs;// - tree[j][nn].L2P_f * L2M_f_QR.solve(tree[j][k].local_rhs);
		}
		else {
			tree[j][nn].multipole_rhs = tree[j][nn].multipole_rhs;// - tree[j][nn].L2P_f * L2M_f_QR.solve(tree[j][k].local_rhs);
		}
	}

	void filings_in_equation_P2P_due_to_neighbors(int j, int k, int nn) {
		// nn: equation number
		// fill-ins between neighbors of k
		for (size_t p = 0; p < 8; p++) {
			int pp = tree[j][k].neighborNumbers[p];//P2P(j,nn)
			if (!tree[j][nn].Eliminated) {
				if(pp != -1) {
					if (!tree[j][pp].Eliminated) { //case 4
						tree[j][nn].P2P[pp] = tree[j][nn].P2P[pp];// - tree[j][nn].L2P_f * L2M_f_QR.solve(tree[j][k].P2M_f[pp]);
						if(is_well_separated(j,nn,pp)) {
							if (nn < pp) { // using symmetry
								if(tree[j][nn].P2P[pp].norm() > 1e-10) {
									// std::cout << "compress_P2P here" << std::endl;
									// compress_P2P(j,nn,pp);
									std::pair<int, int> g(nn,pp);
									P2P.push_back(g);
									// std::cout << "there" << std::endl;
									// std::cout << "j: " << j << "	nn: " << nn << "	pp: " << pp << "	compress_P2P" << std::endl;
								}
							}
							// tree[j][nn].P2P[pp] = Mat(0,0);
						}
					}
					else { //case 2/3
						tree[j][nn].M2P[pp] = tree[j][nn].M2P[pp];// - tree[j][nn].L2P_f * L2M_f_QR.solve(tree[j][k].M2M_f[pp]);
						if(is_well_separated(j,nn,pp)) {
							if (tree[j][nn].M2P[pp].norm() > 1e-10) {
								// std::cout << "compress_M2P here" << std::endl;
								// compress_M2P(j,nn,pp);
								// std::cout << "there" << std::endl;
							}
							// tree[j][nn].M2P[pp] = Mat(0,0);
						}
					}
				}
			}
			else {
				if(pp != -1) {
					if (!tree[j][pp].Eliminated) { //case 2/3
						tree[j][nn].P2L[pp] = tree[j][nn].P2L[pp];// - tree[j][nn].L2P_f * L2M_f_QR.solve(tree[j][k].P2M_f[pp]);
						if(is_well_separated(j,nn,pp)) {
							if (tree[j][nn].P2L[pp].norm() > 1e-10) {
								// std::cout << "compress_P2L here" << std::endl;
								// compress_P2L(j,nn,pp);
								std::pair<int, int> g(nn,pp);
								P2L_M2P.push_back(g);
								// std::cout << "there" << std::endl;
								// std::cout << "j: " << j << "	nn: " << nn << "	pp: " << pp << "	compress_P2L" << std::endl;
							}
							// tree[j][nn].P2L[pp] = Mat(0,0);
						}
					}
					else { // case 1
						tree[j][nn].M2L[pp] = tree[j][nn].M2L[pp];// - tree[j][nn].L2P_f * L2M_f_QR.solve(tree[j][k].M2M_f[pp]);
					}
				}
			}
		}
	}

	void filings_in_equation_P2P_due_to_y(int j, int k, int nn) {
		// fill-ins between k and its neighbors
		if (!tree[j][nn].Eliminated) { //case 2/3
			if(tree[j][nn].M2P[k].size() == 0) {
				tree[j][nn].M2P[k] = Mat(tree[j][nn].L2P_f[k].rows(),tree[j][k].NumMultipoles);//- tree[j][nn].L2P_f * L2M_f_QR.solve(-Mat::Identity(tree[j][k].NumMultipoles, tree[j][k].NumMultipoles));
			}
			else {
				tree[j][nn].M2P[k] = tree[j][nn].M2P[k];// - tree[j][nn].L2P_f * L2M_f_QR.solve(-Mat::Identity(tree[j][k].NumMultipoles, tree[j][k].NumMultipoles));
			}
			if(is_well_separated(j,nn,k)) {
				if (tree[j][nn].M2P[k].norm() > 1e-10) {
					// std::cout << "compress_M2P here" << std::endl;
					// compress_M2P(j,nn,k);
					// std::cout << "there" << std::endl;
				}
				// tree[j][nn].M2P[k] = Mat(0,0);
			}
		}
		else { //case 1
			if (tree[j][nn].M2L[k].size() == 0) {
				tree[j][nn].M2L[k] = Mat(tree[j][nn].L2P_f[k].rows(), tree[j][k].NumMultipoles);//- tree[j][nn].L2P_f * L2M_f_QR.solve(-Mat::Identity(tree[j][k].NumMultipoles, tree[j][k].NumMultipoles));
			}
			else {
				tree[j][nn].M2L[k] = tree[j][nn].M2L[k];// - tree[j][nn].L2P_f * L2M_f_QR.solve(-Mat::Identity(tree[j][k].NumMultipoles, tree[j][k].NumMultipoles));
			}
		}
	}

	///////////////////////////////////////////////
  void rhs_filings_in_equation_M2L_due_to_z(int j, int k, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		rhs_update_in_equation_M2L(j,k,L2M_f_QR);
		// std::cout << "rhs_update_in_equation_M2L done" << std::endl;
	}

	void filings_in_equation_M2L_due_to_z(int j, int k, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		filings_in_equation_M2L_due_to_y(j,k,L2M_f_QR);
		// std::cout << "filings_in_equation_M2L_due_to_y done" << std::endl;
		filings_in_equation_M2L_due_to_neighbors(j,k,L2M_f_QR);
		// std::cout << "filings_in_equation_M2L_due_to_neighbors done" << std::endl;
		// rhs_update_in_equation_M2L(j,k,L2M_f_QR);
		// std::cout << "rhs_update_in_equation_M2L done" << std::endl;
	}

	void rhs_update_in_equation_M2L(int j, int k, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		tree[j][k].multipole_rhs = L2M_f_QR.solve(tree[j][k].local_rhs);
	}

	void filings_in_equation_M2L_due_to_y(int j, int k, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		// std::cout << "tree[j][k].L2M_f.S: " << tree[j][k].L2M_f.rows() << ", " << tree[j][k].L2M_f.cols() << std::endl;
		// std::cout << "tree[j][k].NumMultipoles: " << tree[j][k].NumMultipoles << std::endl;
		tree[j][k].M2L[k] = L2M_f_QR.solve(-Mat::Identity(tree[j][k].NumMultipoles, tree[j][k].NumMultipoles));
		// std::cout << "done" << std::endl;
	}

	void filings_in_equation_M2L_due_to_neighbors(int j, int k, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		for (size_t n = 0; n < 8; n++) {
			int nn = tree[j][k].neighborNumbers[n];//P2P(j,nn)
			if(nn != -1) {
				if (!tree[j][nn].Eliminated) {
					tree[j][k].P2L[nn] = L2M_f_QR.solve(tree[j][k].P2M_f[nn]);
					if(is_well_separated(j,k,nn)) {
						if (tree[j][k].P2L[nn].norm() > 1e-10) {
							// std::cout << "compress_P2L here" << std::endl;
							std::pair<int, int> g(k,nn);
							P2L_M2P.push_back(g);
							// compress_P2L(j,k,nn);
							// std::cout << "there" << std::endl;
							// std::cout << "j: " << j << "	k: " << k << "	nn: " << nn << "	compress_P2L" << std::endl;
						}
						// tree[j][k].P2L[nn] = Mat(0,0);
					}
				}
				else {
					tree[j][k].M2L[nn] = L2M_f_QR.solve(tree[j][k].M2M_f[nn]);
					// if (j==3 && k==6 && nn==10) {
					// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree[j][k].M2L[nn]: " << std::endl << tree[j][k].M2L[nn] << std::endl << std::endl;
					// }
				}
			}
		}
	}

  void rhs_filings_in_equation_P2P_due_to_z(int j, int k, int nn, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
    // std::cout << "rhs_update_in_equation_P2P_due_to_z START" << std::endl;
		rhs_update_in_equation_P2P_due_to_z(j,k,nn,L2M_f_QR);
		// std::cout << "rhs_update_in_equation_P2P_due_to_z done" << std::endl;
	}

	void filings_in_equation_P2P_due_to_z(int j, int k, int nn, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		filings_in_equation_P2P_due_to_y(j,k,nn,L2M_f_QR);
		// std::cout << "filings_in_equation_P2P_due_to_y done" << std::endl;
		filings_in_equation_P2P_due_to_neighbors(j,k,nn,L2M_f_QR);
    // if (k==1 && nn==4) {
    //   std::cout << "tree[j][nn].Eliminated: " << tree[j][nn].Eliminated << std::endl;
    //   std::cout << "tree[j][k].L2M_f: " << tree[j][k].L2M_f.rows() << "," << tree[j][k].L2M_f.cols() << std::endl;
    //   std::cout << "tree[j][nn].L2P_f[k]: " << tree[j][nn].L2P_f[k].rows() << "," << tree[j][nn].L2P_f[k].cols() << std::endl;
    //   std::cout << "tree[j][nn].local_rhs: " << tree[j][k].local_rhs.size() << std::endl;
    //   std::cout << "tree[j][nn].multipole_rhs: " << tree[j][k].multipole_rhs.size() << std::endl;
    //   std::cout << "tree[j][nn].particle_rhs: " << tree[j][k].particle_rhs.size() << std::endl;
    // }
		// std::cout << "filings_in_equation_P2P_due_to_neighbors done" << std::endl;
		// rhs_update_in_equation_P2P_due_to_z(j,k,nn,L2M_f_QR);
		// std::cout << "rhs_update_in_equation_P2P_due_to_z done" << std::endl;
	}

	void rhs_update_in_equation_P2P_due_to_z(int j, int k, int nn, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
    // std::cout << "tree[j][nn].Eliminated: " << tree[j][nn].Eliminated << std::endl;
    // std::cout << "tree[j][k].L2M_f: " << tree[j][k].L2M_f.rows() << "," << tree[j][k].L2M_f.cols() << std::endl;
    // std::cout << "tree[j][nn].L2P_f[k]: " << tree[j][nn].L2P_f[k].rows() << "," << tree[j][nn].L2P_f[k].cols() << std::endl;
    // std::cout << "tree[j][nn].local_rhs: " << tree[j][k].local_rhs.size() << std::endl;
    // std::cout << "tree[j][nn].multipole_rhs: " << tree[j][k].multipole_rhs.size() << std::endl;
    // std::cout << "tree[j][nn].particle_rhs: " << tree[j][k].particle_rhs.size() << std::endl;
		if (!tree[j][nn].Eliminated) {
			tree[j][nn].particle_rhs = tree[j][nn].particle_rhs - tree[j][nn].L2P_f[k] * L2M_f_QR.solve(tree[j][k].local_rhs);
		}
		else {
			tree[j][nn].multipole_rhs = tree[j][nn].multipole_rhs - tree[j][nn].L2P_f[k] * L2M_f_QR.solve(tree[j][k].local_rhs);
		}
	}

	void filings_in_equation_P2P_due_to_neighbors(int j, int k, int nn, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		// nn: equation number
		// fill-ins between neighbors of k
		for (size_t p = 0; p < 8; p++) {
			int pp = tree[j][k].neighborNumbers[p];//P2P(j,nn)
			if (!tree[j][nn].Eliminated) {
				if(pp != -1) {
					if (!tree[j][pp].Eliminated) { //case 4
						tree[j][nn].P2P[pp] = tree[j][nn].P2P[pp] - tree[j][nn].L2P_f[k] * L2M_f_QR.solve(tree[j][k].P2M_f[pp]);
						if(is_well_separated(j,nn,pp)) {
							if (nn < pp) { // using symmetry
								if(tree[j][nn].P2P[pp].norm() > 1e-10) {
									// std::cout << "compress_P2P here" << std::endl;
									// compress_P2P(j,nn,pp);
									std::pair<int, int> g(nn,pp);
									P2P.push_back(g);
									// std::cout << "there" << std::endl;
									// std::cout << "j: " << j << "	nn: " << nn << "	pp: " << pp << "	compress_P2P" << std::endl;
								}
							}
							// tree[j][nn].P2P[pp] = Mat(0,0);
						}
					}
					else { //case 2/3
						tree[j][nn].M2P[pp] = tree[j][nn].M2P[pp] - tree[j][nn].L2P_f[k] * L2M_f_QR.solve(tree[j][k].M2M_f[pp]);
						if(is_well_separated(j,nn,pp)) {
							if (tree[j][nn].M2P[pp].norm() > 1e-10) {
								// std::cout << "compress_M2P here" << std::endl;
								// compress_M2P(j,nn,pp);
								// std::cout << "there" << std::endl;
							}
							// tree[j][nn].M2P[pp] = Mat(0,0);
						}
					}
				}
			}
			else {
				if(pp != -1) {
					if (!tree[j][pp].Eliminated) { //case 2/3
						tree[j][nn].P2L[pp] = tree[j][nn].P2L[pp] - tree[j][nn].L2P_f[k] * L2M_f_QR.solve(tree[j][k].P2M_f[pp]);
						if(is_well_separated(j,nn,pp)) {
							if (tree[j][nn].P2L[pp].norm() > 1e-10) {
								// std::cout << "compress_P2L here" << std::endl;
								// compress_P2L(j,nn,pp);
								std::pair<int, int> g(nn,pp);
								P2L_M2P.push_back(g);
								// std::cout << "there" << std::endl;
								// std::cout << "j: " << j << "	nn: " << nn << "	pp: " << pp << "	compress_P2L" << std::endl;
							}
							// tree[j][nn].P2L[pp] = Mat(0,0);
						}
					}
					else { // case 1
						tree[j][nn].M2L[pp] = tree[j][nn].M2L[pp] - tree[j][nn].L2P_f[k] * L2M_f_QR.solve(tree[j][k].M2M_f[pp]);
						// if (j==3 && nn==6 && pp==10) {
						// 	std::cout << "eliminate_z: k: " << k << std::endl;
						// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree[j][nn].M2L[pp]: " << std::endl << tree[j][nn].M2L[pp] << std::endl << std::endl;
						// }
					}
				}
			}
		}
	}

	void filings_in_equation_P2P_due_to_y(int j, int k, int nn, Eigen::PartialPivLU<Mat>& L2M_f_QR) {
		// fill-ins between k and its neighbors
		if (!tree[j][nn].Eliminated) { //case 2/3
			if(tree[j][nn].M2P[k].size() == 0) {
				tree[j][nn].M2P[k] = - tree[j][nn].L2P_f[k] * L2M_f_QR.solve(-Mat::Identity(tree[j][k].NumMultipoles, tree[j][k].NumMultipoles));
			}
			else {
				tree[j][nn].M2P[k] = tree[j][nn].M2P[k] - tree[j][nn].L2P_f[k] * L2M_f_QR.solve(-Mat::Identity(tree[j][k].NumMultipoles, tree[j][k].NumMultipoles));
			}
			if(is_well_separated(j,nn,k)) {
				if (tree[j][nn].M2P[k].norm() > 1e-10) {
					// std::cout << "compress_M2P here" << std::endl;
					// compress_M2P(j,nn,k);
					// std::cout << "there" << std::endl;
				}
				// tree[j][nn].M2P[k] = Mat(0,0);
			}
		}
		else { //case 1
			if (tree[j][nn].M2L[k].size() == 0) {
				tree[j][nn].M2L[k] = - tree[j][nn].L2P_f[k] * L2M_f_QR.solve(-Mat::Identity(tree[j][k].NumMultipoles, tree[j][k].NumMultipoles));
			}
			else {
				tree[j][nn].M2L[k] = tree[j][nn].M2L[k] - tree[j][nn].L2P_f[k] * L2M_f_QR.solve(-Mat::Identity(tree[j][k].NumMultipoles, tree[j][k].NumMultipoles));
			}
			// if (j==3 && nn==6 && k==10) {
			// 	std::cout << std::endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!! tree[j][nn].M2L[k]: " << std::endl << tree[j][nn].M2L[k] << std::endl << std::endl;
			// }
		}
	}

	void eliminate_phase() {
		// make_basis_unitary();
		// update_M2L_to_unitary_basis();
		for(int j=nLevels; j>=2; --j) {
			if(j != nLevels) {
				initialise_P2P_NonLeafLevel(j);
			}
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				// std::cout << "j: " << j << "	k: " << k  << "------------------------------------" << std::endl;
				// std::cout << "L2P.S: " << tree[j][k].L2P.rows() << ", " << tree[j][k].L2P.cols() << "	P2M.S: " << tree[j][k].P2M.rows() << ", " << tree[j][k].P2M.cols() << std::endl;
				tree[j][k].Eliminated = true;
				eliminate_cluster(j,k);
				for (size_t i = 0; i < P2P.size(); i++) {
					int nn = P2P[i].first;
					int pp = P2P[i].second;
					// std::cout << "compress_P2P	" << nn << ", " << pp << std::endl;
					compress_P2P_Two_Ways(j,nn,pp);
					// std::cout << "DONE" << std::endl;
					tree[j][nn].P2P[pp] = Mat(0,0);
					tree[j][pp].P2P[nn] = Mat(0,0);
				}
				P2P.clear();
				for (size_t i = 0; i < P2L_M2P.size(); i++) {
					int nn = P2L_M2P[i].first;
					int pp = P2L_M2P[i].second;
					// std::cout << "compress_P2L" << pp << ", " << nn << std::endl;
					compress_M2P_P2L(j, pp, nn);
					// std::cout << "DONE" << std::endl;
					tree[j][nn].P2L[pp] = Mat(0,0);
					tree[j][pp].M2P[nn] = Mat(0,0);
				}
				P2L_M2P.clear();
			}
		}
		//solve for particles at level 0; which are nothing but the multipoles at level 1
		int j = 2;//base level
		solve_particles_at_base_level(j);//now multipoles at level 1 are available
	}

	void eliminate_phase_efficient() {
		// make_basis_unitary();
		// update_M2L_to_unitary_basis();
		for(int j=nLevels; j>=2; --j) {
			if(j != nLevels) {
				initialise_P2P_NonLeafLevel(j);
        // std::cout << "-----------DONE-------------------------" << std::endl;
			}
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				// std::cout << "j: " << j << "	k: " << k  << "------------------------------------" << std::endl;

				for (size_t i = 0; i < P2P.size(); i++) {
					int nn = P2P[i].first;
					int pp = P2P[i].second;
					if (k == nn || k == pp) {
						P2P.erase(P2P.begin()+i);
						i--;
						if (tree[j][nn].P2P[pp].size() != 0) {
							// std::cout << "compress_P2P:	" << nn << ", " << pp << std::endl;
							compress_P2P_Two_Ways(j,nn,pp);
							// std::cout << "DONE" << std::endl;
							tree[j][nn].P2P[pp] = Mat(0,0);
							tree[j][pp].P2P[nn] = Mat(0,0);
						}
					}
				}
				// std::cout << "P2P compression DONE" << std::endl;
				for (size_t i = 0; i < P2L_M2P.size(); i++) {
					int nn = P2L_M2P[i].first;
					int pp = P2L_M2P[i].second;
					if (k == nn || k == pp) {
						P2L_M2P.erase(P2L_M2P.begin()+i);
						i--;
						if (tree[j][nn].P2L[pp].size() != 0) {
							// std::cout << "compress_P2L:	" << pp << ", " << nn << std::endl;
							compress_M2P_P2L(j, pp, nn);
							// std::cout << "DONE" << std::endl;
							tree[j][nn].P2L[pp] = Mat(0,0);
							tree[j][pp].M2P[nn] = Mat(0,0);
						}
					}
				}

				// std::cout << "P2L_M2P compression DONE" << std::endl;
				// std::cout << "P2P.size(): " << P2P.size() << std::endl;
				// std::cout << "P2L_M2P.size(): " << P2L_M2P.size() << std::endl;
				// P2P.clear();
				// P2L_M2P.clear();
				tree[j][k].Eliminated = true;
				eliminate_cluster(j,k);
			}
		}
	}

  void rhs_eliminate_phase_efficient() {
    for(int j=nLevels; j>=2; --j) {
      // std::cout << "here" << std::endl;
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
        tree[j][k].Eliminated = false;
      }
      // std::cout << "here" << std::endl;
    }
		for(int j=nLevels; j>=2; --j) {
      if(j != nLevels) {
        assign_particle_rhs_NonLeafLevel(j);
        // std::cout << "-----------DONE-------------------------" << std::endl;
			}
      for (int k=0; k<nBoxesPerLevel[j]; ++k) {
        // std::cout << "k: " << k << std::endl;
        tree[j][k].Eliminated = true;
        rhs_eliminate_cluster(j,k);
			}
		}
	}

	void solve_particles_at_base_level(int j) {
		int NumMultipoles = 0;
		int NumLocals = 0;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			NumMultipoles += tree[j][k].NumMultipoles;
			NumLocals += tree[j][k].NumLocals;
		}
    if (NumLocals == 0 || NumMultipoles == 0) {
      return;
    }
    else {
  		Mat A = Mat::Zero(NumLocals, NumMultipoles);
  		Vec b = Vec::Zero(NumLocals);
  		int r_index = 0;
  		for (int r=0; r<nBoxesPerLevel[j]; ++r) {
  			int s_index = 0;
  			for (int s=0; s<nBoxesPerLevel[j]; ++s) {
  				A.block(r_index, s_index, tree[j][r].NumLocals, tree[j][s].NumMultipoles) = tree[j][r].M2L[s];
  				s_index += tree[j][s].NumMultipoles;
  			}
  			b.segment(r_index, tree[j][r].NumLocals) = tree[j][r].multipole_rhs;
  			r_index += tree[j][r].NumLocals;
  		}
  		// Eigen::PartialPivLU<Mat> A_lu = A.lu();
  		// std::cout << "A.s: " << A.rows() << ", " << A.cols() << std::endl;
  		Eigen::PartialPivLU<Mat> A_qr = A.lu();
  		Vec multipoles = A_qr.solve(b);
  		// std::cout << "A.n: " << A.norm()  << std::endl;
  		// std::cout << "b.n: " << b.norm()  << std::endl;
  		// std::cout << "multipoles.N: " << multipoles.norm() << std::endl;
  		// std::cout << "symmetry err: " << (A-A.transpose()).norm() << std::endl;
  		int index = 0;
  		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
  			tree[j][k].multipoles = multipoles.segment(index, tree[j][k].NumMultipoles);
  			index += tree[j][k].NumMultipoles;
  		}
    }
	}

	void back_substitution_phase() {
    //solve for particles at level 0; which are nothing but the multipoles at level 1
    int j = 2;//base level
    solve_particles_at_base_level(j);//now multipoles at level 1 are available
		for (size_t j = 2; j <= nLevels; j++) {
			for (int k=nBoxesPerLevel[j]-1; k>=0; k--) {
				// obtain x and z using P2P and U.transpose() equations
				// Eigen::PartialPivLU<Mat> P2P_lu = tree[j][k].P2P[k].lu();
				// std::cout << "j: " << j << "	k: " << k << "	------------------------" << std::endl;
				// std::cout << "tree[j][k].P2P[k].S: " << tree[j][k].P2P[k].rows() << ", " << tree[j][k].P2P[k].cols() << std::endl;
				if (tree[j][k].P2P[k].size() != 0) {
					Eigen::PartialPivLU<Mat> P2P_qr = tree[j][k].P2P[k].lu();
					Vec rhs1;
					rhs1 = tree[j][k].particle_rhs;
					for (size_t n = 0; n < 8; n++) {
						int nn = tree[j][k].neighborNumbers[n];
						if (nn != -1) {
							if (!tree[j][nn].Eliminated) {
								rhs1 = rhs1 - tree[j][k].P2P[nn]*tree[j][nn].particles;
							}
							else {
								rhs1 = rhs1 - tree[j][k].M2P[nn]*tree[j][nn].multipoles;
							}
						}
					}
					Vec rhs2 = tree[j][k].multipoles - tree[j][k].P2M * (P2P_qr .solve(rhs1));
					Mat Atemp = -tree[j][k].P2M * (P2P_qr .solve(tree[j][k].L2P));
					// std::cout << "Atemp.S: " << Atemp.rows() << ", " << Atemp.cols() << std::endl;
					if (Atemp.size() != 0) {
						Eigen::PartialPivLU<Mat> Atemp_qr = Atemp.lu();
						tree[j][k].locals = Atemp_qr .solve(rhs2);
					}
					else {
						tree[j][k].locals = Vec(0);
					}
					// tree[j][k].locals = Atemp.lu() .solve(rhs2);
					tree[j][k].particles = P2P_qr .solve(rhs1-tree[j][k].L2P * tree[j][k].locals);
				}
				else {
					tree[j][k].particles = Vec(0);
				}
				tree[j][k].Eliminated = false;
				if (nLevels != j) {
					int index = 0;
					for (size_t c = 0; c < 4; c++) { // get y^{k+1} from x^{k}
						tree[j+1][4*k+c].multipoles = tree[j][k].particles.segment(index, tree[j+1][4*k+c].NumMultipoles);
						index += tree[j+1][4*k+c].NumMultipoles;
					}
				}
			}
		}
	}

	void check2() {
		std::cout << "multipoles and locals of leaf level boxes" << std::endl << std::endl;
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			std::cout << "k: " << k << ";	" << tree[nLevels][k].NumMultipoles << ", " << tree[nLevels][k].NumLocals << ", " << tree[nLevels][k].NumMultipoles-tree[nLevels][k].NumLocals << std::endl;
		}
	}

  void getx(Vec &x) {
    int indexVec[nBoxesPerLevel[nLevels]];
    for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
      if(k==0) {
        indexVec[k] = 0;
      }
      else {
        indexVec[k] = indexVec[k-1] + tree[nLevels][k-1].chargeLocations.size();
      }
    }
    x = Vec::Zero(gridPoints.size()); //all particles
    // #pragma omp parallel for
    for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
      x.segment(indexVec[k], tree[nLevels][k].particles.size()) = tree[nLevels][k].particles;
      // index += tree[nLevels][k].particles.size();
    }
  }

  double error_check() {
    // true matrix A generation
/////////////////////////////////
    // Mat A(N,N);
    // int index1 = 0;
    // int index2 = 0;
    // for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    // 	for (size_t b = 0; b < nBoxesPerLevel[nLevels]; b++) {
    // 		A.block(index1, index2, tree[nLevels][k].chargeLocations.size(), tree[nLevels][b].chargeLocations.size()) = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][b].chargeLocations);
    // 		index2 += tree[nLevels][b].chargeLocations.size();
    // 	}
    // 	index2 = 0;
    // 	index1 += tree[nLevels][k].chargeLocations.size();
    // }
    /////////////////////////////////
    int indexVec[nBoxesPerLevel[nLevels]];
    for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
      if(k==0) {
        indexVec[k] = 0;
      }
      else {
        indexVec[k] = indexVec[k-1] + tree[nLevels][k-1].chargeLocations.size();
      }
    }

    Vec x = Vec::Zero(gridPoints.size()); //all particles
    // int index = 0;
    #pragma omp parallel for
    for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
      x.segment(indexVec[k], tree[nLevels][k].particles.size()) = tree[nLevels][k].particles;
      // index += tree[nLevels][k].particles.size();
    }
    // index = 0;
    b_error_check = Vec::Zero(N);
    #pragma omp parallel for
    for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
      b_error_check.segment(indexVec[k], tree[nLevels][k].chargeLocations.size()) = tree[nLevels][k].particle_rhs;
      // index += tree[nLevels][k].chargeLocations.size();
    }
    // the following mimics ''Vec rhs = A*x;''; done to avoid memory issue

    Vec rhs = Vec::Zero(gridPoints.size()); //all particles
    /////////////////////////////////
    // Vec rhs3 = Vec::Zero(gridPoints.size()); //all particles
    // index = 0;
    // for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    //   int index2 = 0;
    //   for (size_t b = 0; b < nBoxesPerLevel[nLevels]; b++) {
    //     Mat B = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][b].chargeLocations);
    //     rhs3.segment(index, tree[nLevels][k].particles.size()) += B*tree[nLevels][b].particles;
    //     index2 += tree[nLevels][b].chargeLocations.size();
    //   }
    //   index += tree[nLevels][k].particles.size();
    // }
    /////////////////////////////////

    #pragma omp parallel for
    for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
      // #pragma omp parallel for
      for (size_t b = 0; b < nBoxesPerLevel[nLevels]; b++) {
        Mat B = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][b].chargeLocations);
        rhs.segment(indexVec[k], tree[nLevels][k].chargeLocations.size()) += B*tree[nLevels][b].particles;
      }
    }

    // Vec rhs2 = A*x;

    // apprximate x to be saved: x
    // RHS to be saved: Ax

    // string fname;
    // fname = "rhs_x/rhs_real_" + std::to_string((int)Qchoice) + "_" + std::to_string((int)N) + ".bin";
    // storedata::save_vec(fname, rhs.real());
    // fname = "rhs_x/rhs_imag_" + std::to_string((int)Qchoice) + "_" + std::to_string((int)N) + ".bin";
    // storedata::save_vec(fname, rhs.imag());
    // fname = "rhs_x/x_real_" + std::to_string((int)Qchoice) + "_" + std::to_string((int)N) + ".bin";
    // storedata::save_vec(fname, x.real());
    // fname = "rhs_x/x_imag_" + std::to_string((int)Qchoice) + "_" + std::to_string((int)N) + ".bin";
    // storedata::save_vec(fname, x.imag());

    // string fname;
    // fname = "rhs_x/rhs_real_" + std::to_string((int)Qchoice) + "_" + std::to_string((int)kappa) + "_" + std::to_string((int)N) + ".bin";
    // storedata::save_vec(fname, rhs.real());
    // fname = "rhs_x/rhs_imag_" + std::to_string((int)Qchoice) + "_" + std::to_string((int)kappa) + "_" + std::to_string((int)N) + ".bin";
    // storedata::save_vec(fname, rhs.imag());
    // fname = "rhs_x/x_real_" + std::to_string((int)Qchoice) + "_" + std::to_string((int)kappa) + "_" + std::to_string((int)N) + ".bin";
    // storedata::save_vec(fname, x.real());
    // fname = "rhs_x/x_imag_" + std::to_string((int)Qchoice) + "_" + std::to_string((int)kappa) + "_" + std::to_string((int)N) + ".bin";
    // storedata::save_vec(fname, x.imag());

    Vec E = rhs-b_error_check;
    std::cout << "gen rhs x" << std::endl;
    std::cout << "abs error: " << E.lpNorm<Infinity>() << std::endl;
    std::cout << "b norm: " << b_error_check.lpNorm<Infinity>() << std::endl;
    return E.norm()/b_error_check.norm();
  }

	void outputCropGraphExtendedSparseMatrixLevel2(std::string filename) {
		std::ofstream myfile;
		myfile.open(filename.c_str());
		myfile << "\\documentclass{article}" << std::endl;
		myfile << "\\usepackage{geometry}" << std::endl;
 		myfile << "\\geometry{left=20mm, top=20mm,papersize={120cm,58cm},}" << std::endl;
		myfile << "\\usepackage{pgfplots}" << std::endl;
		myfile << "\\usepackage{amsmath}" << std::endl;
		myfile << "\\usepackage{amssymb}" << std::endl;
		myfile << "\\usetikzlibrary{arrows}" << std::endl;
		myfile << "\\usetikzlibrary{shapes}" << std::endl;
		myfile << "\\usetikzlibrary{arrows.meta,calc,positioning,quotes}" << std::endl;
		myfile << "\\pgfplotsset{compat=newest,width=300cm, height=50cm,xtick pos=left, ytick pos=left,}" << std::endl;
		myfile << "\\begin{document}" << std::endl;
		myfile << "\\begin{tikzpicture}[->,>=stealth',auto,node distance=5cm, thick,main node/.style={circle,draw,minimum size=1cm,font=\\sffamily\\Large\\bfseries}]" << std::endl;
		myfile << "\\node[main node] (0)  {\\scalebox{2.5}{$x_{0}^{(3)}$}};" << std::endl;
		// myfile << "\\node[main node] (0)  {0};" << std::endl;
		for (size_t i = 1; i < nBoxesPerLevel[nLevels]; i++) {
			if(i < 4)
				myfile << "\\node[main node] (";
			else
				myfile << "\\node[main node,color=white] (";
			myfile << i;
			myfile << ")   [right of=";
			myfile << i-1;
			myfile << "]  {\\scalebox{2.5}{$\\displaystyle x_{";
			myfile << i;
			myfile << "}^{(3)}$}};" << std::endl;
		}
		//self
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < 4; i++) {
			myfile << "(";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=red, loop below] node [color=red, pos=0.5, above] {} (";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;

		//P2P
		for (size_t i = 0; i < 4; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for (size_t n = 0; n < 8; n++) {//neighbors of k
				int nn = tree[nLevels][i].neighborNumbers[n];//P2P(j,nn)
				if (nn != -1) {
					if (i < nn) {//avoiding two sided
						myfile << "(";
						myfile << i;
						if (i%2 == 0) {
							myfile << ") edge[line width=2pt, draw=red, bend right, <->] node [color=red, pos=0.5, above] {} (";
						}
						else {
							myfile << ") edge[line width=2pt, draw=red, bend left, <->] node [color=red, pos=0.5, above] {} (";
						}
						myfile << nn;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		//z3
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			if(i<4)
				myfile << "\\node[main node] (z";
			else
				myfile << "\\node[main node,color=white] (z";
			myfile << i;
			myfile << ") [above of=";
			myfile << i;
			myfile << ",yshift=2cm] {\\scalebox{2.5}{$\\displaystyle z_{";
			myfile << i;
			myfile << "}^{(3)}$}};" << std::endl;
		}
		//x3<->z3
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < 4; i++) {
			myfile << "(";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=green, <->] node [color=green, pos=0.5, right] {} (z";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;
		//y3
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			if(i<4)
				myfile << "\\node[main node] (y";
			else
				myfile << "\\node[main node,color=white] (y";
			myfile << i;
			myfile << ") [above of=z";
			myfile << i;
			// myfile << ",yshift=6cm] {};";
			myfile << ",yshift=2.5cm] {\\scalebox{2.5}{$\\displaystyle y_{";
			myfile << i;
			myfile << "}^{(3)}$}};" << std::endl;
		}
		//z3<->y3
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < 4; i++) {
			myfile << "(z";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=blue, <->] node [color=green, pos=0.5, right] {} (y";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;
		//M2L3
		for (size_t i = 0; i < 4; i++) {
		// for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<16; ++in) {
				if(tree[nLevels][i].innerNumbers[in] != -1) {
					int kIL = tree[nLevels][i].innerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(y";
						myfile << i;
						if (kIL%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		for (size_t i = 0; i < 4; i++) {
		// for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<24; ++in) {
				if(tree[nLevels][i].outerNumbers[in] != -1) {
					int kIL = tree[nLevels][i].outerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(y";
						myfile << i;
						if (kIL%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		// //z2
		// int u=0;
		// for (size_t i = 1; i < nBoxesPerLevel[nLevels]; i+=4) {
		// 	// \node[main node, above=of $(8.north)!0.5!(9.north)$] (10) {10};
		// 	if (i<4)
		// 		myfile << "\\node[main node, above=of $(y";
		// 	else
		// 		myfile << "\\node[main node, color=white, above=of $(y";
		// 	myfile << i;
		// 	myfile << ".north)!0.5!(y";
		// 	myfile << i+1;
		// 	myfile << ".north)$,yshift=2cm] (2z";
		// 	myfile << u;
		// 	myfile << ") {\\scalebox{2.5}{$\\displaystyle z_{";
		// 	myfile << u;
		// 	myfile << "}^{(2)}$}};" << std::endl;
		// 	u++;
		// }
		// // //y2<->z2
		// myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		// for (size_t i = 0; i < 4; i++) {
		// 	myfile << "(y";
		// 	myfile << i;
		// 	myfile << ") edge[line width=2pt, draw=green, <->] node [color=green, pos=0.5, right] {} (2z";
		// 	myfile << int(i/4);
		// 	myfile << ")" << std::endl;
		// }
		// myfile << ";" << std::endl;
		// // //z2<->y2
		// for (size_t i = 0; i < nBoxesPerLevel[nLevels-1]; i++) {
		// 	if(i<1)
		// 		myfile << "\\node[main node] (2y";
		// 	else
		// 		myfile << "\\node[main node, color=white] (2y";
		// 	myfile << i;
		// 	myfile << ") [above of=2z";
		// 	myfile << i;
		// 	// myfile << ",yshift=6cm] {};";
		// 	myfile << ",yshift=2cm] {\\scalebox{2.5}{$\\displaystyle y_{";
		// 	myfile << i;
		// 	myfile << "}^{(2)}$}};" << std::endl;
		// }
		// // //z3<->y3
		// myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		// for (size_t i = 0; i < 1; i++) {
		// 	myfile << "(2z";
		// 	myfile << i;
		// 	myfile << ") edge[line width=2pt, draw=blue, <->] node [color=green, pos=0.5, right] {} (2y";
		// 	myfile << i;
		// 	myfile << ")" << std::endl;
		// }
		// myfile << ";" << std::endl;
		//
		// //M2L2
		// for (size_t i = 0; i < 1; i++) {
		// 	myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		// 	for(int in=0; in<16; ++in) {
		// 		if(tree[nLevels-1][i].innerNumbers[in] != -1) {
		// 			int kIL = tree[nLevels-1][i].innerNumbers[in];
		// 			if (i < kIL) {//avoiding two sided
		// 				myfile << "(2y";
		// 				myfile << i;
		// 				if (kIL%2 == 0) {
		// 					myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (2y";
		// 				}
		// 				else {
		// 					myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (2y";
		// 				}
		// 				myfile << kIL;
		// 				myfile << ")" << std::endl;
		// 			}
		// 		}
		// 	}
		// 	myfile << ";" << std::endl;
		// }
		// for (size_t i = 0; i < 1; i++) {
		// 	myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		// 	for(int in=0; in<24; ++in) {
		// 		if(tree[nLevels-1][i].outerNumbers[in] != -1) {
		// 			int kIL = tree[nLevels-1][i].outerNumbers[in];
		// 			if (i < kIL) {//avoiding two sided
		// 				myfile << "(2y";
		// 				myfile << i;
		// 				if (kIL%2 == 0) {
		// 					myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (2y";
		// 				}
		// 				else {
		// 					myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (2y";
		// 				}
		// 				myfile << kIL;
		// 				myfile << ")" << std::endl;
		// 			}
		// 		}
		// 	}
		// 	myfile << ";" << std::endl;
		// }
		myfile << "\\end{tikzpicture}" << std::endl;
		myfile << "\\end{document}" << std::endl;
		}

		void outputGraphMatrixLevel2(std::string filename) {
			std::ofstream myfile;
			myfile.open(filename.c_str());
			myfile << "\\documentclass{article}" << std::endl;
			myfile << "\\usepackage{geometry}" << std::endl;
	 		myfile << "\\geometry{left=20mm, top=20mm,papersize={120cm,60cm},}" << std::endl;
			myfile << "\\usepackage{pgfplots}" << std::endl;
			myfile << "\\usepackage{amsmath}" << std::endl;
			myfile << "\\usepackage{amssymb}" << std::endl;
			myfile << "\\usetikzlibrary{arrows}" << std::endl;
			myfile << "\\usetikzlibrary{shapes}" << std::endl;
			myfile << "\\usetikzlibrary{arrows.meta,calc,positioning,quotes}" << std::endl;
			myfile << "\\pgfplotsset{compat=newest,width=240cm, height=50cm,xtick pos=left, ytick pos=left,}" << std::endl;
			myfile << "\\begin{document}" << std::endl;
			myfile << "\\begin{tikzpicture}[->,>=stealth',auto,node distance=1cm, thick,main node/.style={circle,draw,minimum size=.1cm,font=\\sffamily\\Large\\bfseries}]" << std::endl;
			myfile << "\\node[main node] (0)  {};" << std::endl;
			// myfile << "\\node[main node] (0)  {0};" << std::endl;
			for (size_t i = 1; i < nBoxesPerLevel[nLevels]; i++) {
				myfile << "\\node[main node] (";
				myfile << i;
				myfile << ")   [right of=";
				myfile << i-1;
				myfile << "]  {";
				// myfile << i;
				myfile << "};" << std::endl;
			}
			//self
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
				myfile << "(";
				myfile << i;
				myfile << ") edge[line width=2pt, draw=red, loop below] node [color=red, pos=0.5, above] {} (";
				myfile << i;
				myfile << ")" << std::endl;
			}
			myfile << ";" << std::endl;

			//P2P
			for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
				myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
				for (size_t n = 0; n < nBoxesPerLevel[nLevels]; n++) {//neighbors of k
					// int nn = tree[nLevels][i].neighborNumbers[n];//P2P(j,nn)
					if (n != i) {
						if (i < n) {//avoiding two sided
							myfile << "(";
							myfile << i;
							if (i%2 == 0) {
								myfile << ") edge[line width=2pt, draw=red, bend right, <->] node [color=red, pos=0.5, above] {} (";
							}
							else {
								myfile << ") edge[line width=2pt, draw=red, bend left, <->] node [color=red, pos=0.5, above] {} (";
							}
							myfile << n;
							myfile << ")" << std::endl;
						}
					}
				}
				myfile << ";" << std::endl;
			}
			myfile << "\\end{tikzpicture}" << std::endl;
			myfile << "\\end{document}" << std::endl;
			}

	void outputGraphExtendedSparseMatrixLevel2(std::string filename) {
		std::ofstream myfile;
		myfile.open(filename.c_str());
		myfile << "\\documentclass{article}" << std::endl;
		myfile << "\\usepackage{geometry}" << std::endl;
 		myfile << "\\geometry{left=20mm, top=20mm,papersize={120cm,60cm},}" << std::endl;
		myfile << "\\usepackage{pgfplots}" << std::endl;
		myfile << "\\usepackage{amsmath}" << std::endl;
		myfile << "\\usepackage{amssymb}" << std::endl;
		myfile << "\\usetikzlibrary{arrows}" << std::endl;
		myfile << "\\usetikzlibrary{shapes}" << std::endl;
		myfile << "\\usetikzlibrary{arrows.meta,calc,positioning,quotes}" << std::endl;
		myfile << "\\pgfplotsset{compat=newest,width=240cm, height=50cm,xtick pos=left, ytick pos=left,}" << std::endl;
		myfile << "\\begin{document}" << std::endl;
		myfile << "\\begin{tikzpicture}[->,>=stealth',auto,node distance=1cm, thick,main node/.style={circle,draw,minimum size=.1cm,font=\\sffamily\\Large\\bfseries}]" << std::endl;
		myfile << "\\node[main node] (0)  {};" << std::endl;
		// myfile << "\\node[main node] (0)  {0};" << std::endl;
		for (size_t i = 1; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "\\node[main node] (";
			myfile << i;
			myfile << ")   [right of=";
			myfile << i-1;
			myfile << "]  {";
			// myfile << i;
			myfile << "};" << std::endl;
		}
		//self
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "(";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=red, loop below] node [color=red, pos=0.5, above] {} (";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;

		//P2P
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for (size_t n = 0; n < 8; n++) {//neighbors of k
				int nn = tree[nLevels][i].neighborNumbers[n];//P2P(j,nn)
				if (nn != -1) {
					if (i < nn) {//avoiding two sided
						myfile << "(";
						myfile << i;
						if (i%2 == 0) {
							myfile << ") edge[line width=2pt, draw=red, bend right, <->] node [color=red, pos=0.5, above] {} (";
						}
						else {
							myfile << ") edge[line width=2pt, draw=red, bend left, <->] node [color=red, pos=0.5, above] {} (";
						}
						myfile << nn;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		//z3
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "\\node[main node] (z";
			myfile << i;
			myfile << ") [above of=";
			myfile << i;
			myfile << ",yshift=6cm] {};" << std::endl;
		}
		//x3<->z3
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "(";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=green, <->] node [color=green, pos=0.5, right] {} (z";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;
		//y3
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "\\node[main node] (y";
			myfile << i;
			myfile << ") [above of=z";
			myfile << i;
			// myfile << ",yshift=6cm] {};";
			myfile << ",yshift=6cm] {};" << std::endl;
		}
		//z3<->y3
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "(z";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=blue, <->] node [color=green, pos=0.5, right] {} (y";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;
		//M2L3
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
		// for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<16; ++in) {
				if(tree[nLevels][i].innerNumbers[in] != -1) {
					int kIL = tree[nLevels][i].innerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(y";
						myfile << i;
						if (i%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
		// for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<24; ++in) {
				if(tree[nLevels][i].outerNumbers[in] != -1) {
					int kIL = tree[nLevels][i].outerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(y";
						myfile << i;
						if (i%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		// //z2
		// int u=0;
		// for (size_t i = 1; i < nBoxesPerLevel[nLevels]; i+=4) {
		// 	// \node[main node, above=of $(8.north)!0.5!(9.north)$] (10) {10};
		// 	myfile << "\\node[main node, above=of $(y";
		// 	myfile << i;
		// 	myfile << ".north)!0.5!(y";
		// 	myfile << i+1;
		// 	myfile << ".north)$,yshift=6cm] (2z";
		// 	myfile << u;
		// 	myfile << ") {};" << std::endl;
		// 	u++;
		// }
		// //y2<->z2
		// myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		// for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
		// 	myfile << "(y";
		// 	myfile << i;
		// 	myfile << ") edge[line width=2pt, draw=green, <->] node [color=green, pos=0.5, right] {} (2z";
		// 	myfile << int(i/4);
		// 	myfile << ")" << std::endl;
		// }
		// myfile << ";" << std::endl;
		// //z2<->y2
		// for (size_t i = 0; i < nBoxesPerLevel[nLevels-1]; i++) {
		// 	myfile << "\\node[main node] (2y";
		// 	myfile << i;
		// 	myfile << ") [above of=2z";
		// 	myfile << i;
		// 	// myfile << ",yshift=6cm] {};";
		// 	myfile << ",yshift=6cm] {};" << std::endl;
		// }
		// //z3<->y3
		// myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		// for (size_t i = 0; i < nBoxesPerLevel[nLevels-1]; i++) {
		// 	myfile << "(2z";
		// 	myfile << i;
		// 	myfile << ") edge[line width=2pt, draw=blue, <->] node [color=green, pos=0.5, right] {} (2y";
		// 	myfile << i;
		// 	myfile << ")" << std::endl;
		// }
		// myfile << ";" << std::endl;
		//
		// //M2L2
		// for (size_t i = 0; i < nBoxesPerLevel[nLevels-1]; i++) {
		// // for (size_t i = 0; i < 1; i++) {
		// 	myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		// 	for(int in=0; in<16; ++in) {
		// 		if(tree[nLevels-1][i].innerNumbers[in] != -1) {
		// 			int kIL = tree[nLevels-1][i].innerNumbers[in];
		// 			if (i < kIL) {//avoiding two sided
		// 				myfile << "(2y";
		// 				myfile << i;
		// 				if (i%2 == 0) {
		// 					myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (2y";
		// 				}
		// 				else {
		// 					myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (2y";
		// 				}
		// 				myfile << kIL;
		// 				myfile << ")" << std::endl;
		// 			}
		// 		}
		// 	}
		// 	myfile << ";" << std::endl;
		// }
		// for (size_t i = 0; i < nBoxesPerLevel[nLevels-1]; i++) {
		// // for (size_t i = 0; i < 1; i++) {
		// 	myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		// 	for(int in=0; in<24; ++in) {
		// 		if(tree[nLevels-1][i].outerNumbers[in] != -1) {
		// 			int kIL = tree[nLevels-1][i].outerNumbers[in];
		// 			if (i < kIL) {//avoiding two sided
		// 				myfile << "(2y";
		// 				myfile << i;
		// 				if (i%2 == 0) {
		// 					myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (2y";
		// 				}
		// 				else {
		// 					myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (2y";
		// 				}
		// 				myfile << kIL;
		// 				myfile << ")" << std::endl;
		// 			}
		// 		}
		// 	}
		// 	myfile << ";" << std::endl;
		// }
		myfile << "\\end{tikzpicture}" << std::endl;
		myfile << "\\end{document}" << std::endl;
		}

	void outputExtendedSparseMatrixLevel2(std::string filename) {
		double xcenter = 0.0;
		double ycenter = 0.0;
		double Lx = L;
		double Ly = L;

		std::ofstream myfile;
		myfile.open(filename.c_str());
		myfile << "\\documentclass{standalone}" << std::endl;
		myfile << "\\usepackage{tikz}" << std::endl;
		myfile << "\\begin{document}" << std::endl;
		myfile << "\\begin{tikzpicture}" << std::endl;
		// x^2
		int j=nLevels;
		int sizeOfMatrix = 3*nBoxesPerLevel[nLevels];
		for (size_t i = nLevels-1; i >= 2; i--) {
			sizeOfMatrix += 2*nBoxesPerLevel[i];
		}
		std::cout << "sizeOfMatrix: " << sizeOfMatrix << std::endl;
		double x,x2;
		double y;
		double diameter = 10.0/sizeOfMatrix;
		// myfile << "\\draw[step=";
		// myfile << diameter;
		// myfile << "] (0,0) grid +(10,10);" << std::endl;
		myfile << "\\draw (0,0) rectangle +(10,10);" << std::endl;

		// x3
		for (int i=0; i<nBoxesPerLevel[nLevels]; ++i) {//rows
			//self
			x = 0.0;
			y = 10.0-diameter;
			for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
				myfile << "\\draw[draw=none,fill=red] (" << x << ",";
				myfile << y << ") rectangle +(";
				myfile << diameter << ",";
				myfile << diameter << ");" << std::endl;
				x += 2*diameter;
				y -= 2*diameter;
			}
			//neighbors
			for (size_t n = 0; n < 8; n++) {//neighbors of k
				int nn = tree[j][i].neighborNumbers[n];//P2P(j,nn)
				if (nn != -1) {
					double diameter = 10.0/sizeOfMatrix;
					double kx = i*2*diameter;
					double ky = 10-(i+1)*2*diameter;
					double nnx = nn*2*diameter;
					double nny = 10-(nn+1)*2*diameter+diameter;
					myfile << "\\draw[draw=none,fill=red] (" << kx << ",";
					myfile << nny << ") rectangle +(";
					myfile << diameter << ",";
					myfile << diameter << ");" << std::endl;
				}
			}
		}
			//L2P
			x = diameter;
			y = 10.0-diameter;
			for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
				myfile << "\\draw[draw=none,fill=green] (" << x << ",";
				myfile << y << ") rectangle +(";
				myfile << diameter << ",";
				myfile << diameter << ");" << std::endl;
				x += 2*diameter;
				y -= 2*diameter;
			}
			//z3
			// for (int i=0; i<nBoxesPerLevel[nLevels]; ++i) {//rows
				//V
				x = 0.0;
				x2 = 2*nBoxesPerLevel[nLevels]*diameter;
				y = 10.0-2*diameter;
				for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
					myfile << "\\draw[draw=none,fill=green] (" << x << ",";
					myfile << y << ") rectangle +(";
					myfile << diameter << ",";
					myfile << diameter << ");" << std::endl;
					x += 2*diameter;
					//-I
					myfile << "\\draw[draw=none,fill=blue] (" << x2 << ",";
					myfile << y << ") rectangle +(";
					myfile << diameter << ",";
					myfile << diameter << ");" << std::endl;
					// if ((k+1)%4 == 0)
					// 	x2 += 2*diameter;
					// else
					x2 += diameter;
					y -= 2*diameter;
				}
			//y3 = x2
			x = diameter;
			y = 10.0-2*diameter*nBoxesPerLevel[nLevels]-diameter;
			for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
				//-I
				myfile << "\\draw[draw=none,fill=blue] (" << x << ",";
				myfile << y << ") rectangle +(";
				myfile << diameter << ",";
				myfile << diameter << ");" << std::endl;
				x += 2*diameter;
				// if ((k+1)%4 == 0)
				// 	y -= 2*diameter;
				// else
				y -= diameter;
			}
			for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
				//IL
				for(int in=0; in<16; ++in) {
					if(tree[j][k].innerNumbers[in] != -1) {
						int kIL = tree[j][k].innerNumbers[in];
						// double kx = kIL*diameter + 2*diameter*nBoxesPerLevel[nLevels] + int(kIL/4)*diameter;
						// double nny = 10-((k+1)*diameter + 2*diameter*nBoxesPerLevel[nLevels]) - int(k/4)*diameter;
						double kx = kIL*diameter + 2*diameter*nBoxesPerLevel[nLevels];
						double nny = 10-((k+1)*diameter + 2*diameter*nBoxesPerLevel[nLevels]);
						myfile << "\\draw[draw=none,fill=orange] (" << kx << ",";
						myfile << nny << ") rectangle +(";
						myfile << diameter << ",";
						myfile << diameter << ");" << std::endl;
					}
				}
				for(int on=0; on<24; ++on) {
					if(tree[j][k].outerNumbers[on] != -1) {
						int kIL = tree[j][k].outerNumbers[on];
						// double kx = kIL*diameter + 2*diameter*nBoxesPerLevel[nLevels] + int(kIL/4)*diameter;
						// double nny = 10-((k+1)*diameter+ 2*diameter*nBoxesPerLevel[nLevels]) - int(k/4)*diameter;
						double kx = kIL*diameter + 2*diameter*nBoxesPerLevel[nLevels];
						double nny = 10-((k+1)*diameter+ 2*diameter*nBoxesPerLevel[nLevels]);
						myfile << "\\draw[draw=none,fill=orange] (" << kx << ",";
						myfile << nny << ") rectangle +(";
						myfile << diameter << ",";
						myfile << diameter << ");" << std::endl;
					}
				}
			}
			// z2
			// x = 2*diameter*nBoxesPerLevel[nLevels];
			// y = 10-(4*diameter+ 2*diameter*nBoxesPerLevel[nLevels])-diameter;
			// x2 = 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1];
			// double y2 = 10-(4*diameter+ 2*diameter*nBoxesPerLevel[nLevels])-diameter;
			// for (int k=0; k<nBoxesPerLevel[nLevels-1]; ++k) {
			// 	myfile << "\\draw[draw=none,fill=green] (" << x << ",";
			// 	myfile << y << ") rectangle +(";
			// 	myfile << 4*diameter << ",";
			// 	myfile << diameter << ");" << std::endl;
			// 	x += 5*diameter;
			// 	y -= 5*diameter;
			// 	//-I
			// 	myfile << "\\draw[draw=none,fill=blue] (" << x2 << ",";
			// 	myfile << y2 << ") rectangle +(";
			// 	myfile << diameter << ",";
			// 	myfile << diameter << ");" << std::endl;
			// 	x2 += diameter;
			// 	y2 -= 5*diameter;
			// }
			// y = 10-(2*diameter*nBoxesPerLevel[nLevels])-4*diameter;
			// x = 2*diameter*nBoxesPerLevel[nLevels]+4*diameter;
			// for (int k=0; k<nBoxesPerLevel[nLevels-1]; ++k) {
			// 	myfile << "\\draw[draw=none,fill=green] (" << x << ",";
			// 	myfile << y << ") rectangle +(";
			// 	myfile << diameter << ",";
			// 	myfile << 4*diameter << ");" << std::endl;
			// 	x += 5*diameter;
			// 	y -= 5*diameter;
			// }
			//y2
			// j = nLevels-1;
			// x = 2*diameter*nBoxesPerLevel[nLevels] + 4*diameter;
			// for (int k=0; k<nBoxesPerLevel[nLevels-1]; ++k) {
			// 	double nny = 10-((k+1)*diameter + 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1]);
			// 	myfile << "\\draw[draw=none,fill=blue] (" << x << ",";
			// 	myfile << nny << ") rectangle +(";
			// 	myfile << diameter << ",";
			// 	myfile << diameter << ");" << std::endl;
			// 	x += 4*diameter;
			// }
			// for (int k=0; k<nBoxesPerLevel[nLevels-1]; ++k) {
			// 	//IL
			// 	for(int in=0; in<16; ++in) {
			// 		if(tree[j][k].innerNumbers[in] != -1) {
			// 			int kIL = tree[j][k].innerNumbers[in];
			// 			double kx = kIL*diameter + 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1];
			// 			double nny = 10-((k+1)*diameter + 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1]);
			// 			myfile << "\\draw[draw=none,fill=orange] (" << kx << ",";
			// 			myfile << nny << ") rectangle +(";
			// 			myfile << diameter << ",";
			// 			myfile << diameter << ");" << std::endl;
			// 		}
			// 	}
			// 	for(int on=0; on<24; ++on) {
			// 		if(tree[j][k].outerNumbers[on] != -1) {
			// 			int kIL = tree[j][k].outerNumbers[on];
			// 			double kx = kIL*diameter + 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1];
			// 			double nny = 10-((k+1)*diameter + 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1]);
			// 			myfile << "\\draw[draw=none,fill=orange] (" << kx << ",";
			// 			myfile << nny << ") rectangle +(";
			// 			myfile << diameter << ",";
			// 			myfile << diameter << ");" << std::endl;
			// 		}
			// 	}
			// }

		myfile << "\\end{tikzpicture}" << std::endl;
		myfile << "\\end{document}" << std::endl;
		myfile.close();
	}

	void outputCropGraphExtendedSparseMatrix(std::string filename) {
		std::ofstream myfile;
		myfile.open(filename.c_str());
		myfile << "\\documentclass{article}" << std::endl;
		myfile << "\\usepackage{geometry}" << std::endl;
 		myfile << "\\geometry{left=20mm, top=20mm,papersize={120cm,58cm},}" << std::endl;
		myfile << "\\usepackage{pgfplots}" << std::endl;
		myfile << "\\usepackage{amsmath}" << std::endl;
		myfile << "\\usepackage{amssymb}" << std::endl;
		myfile << "\\usetikzlibrary{arrows}" << std::endl;
		myfile << "\\usetikzlibrary{shapes}" << std::endl;
		myfile << "\\usetikzlibrary{arrows.meta,calc,positioning,quotes}" << std::endl;
		myfile << "\\pgfplotsset{compat=newest,width=300cm, height=50cm,xtick pos=left, ytick pos=left,}" << std::endl;
		myfile << "\\begin{document}" << std::endl;
		myfile << "\\begin{tikzpicture}[->,>=stealth',auto,node distance=5cm, thick,main node/.style={circle,draw,minimum size=1cm,font=\\sffamily\\Large\\bfseries}]" << std::endl;
		myfile << "\\node[main node] (0)  {\\scalebox{2.5}{$x_{0}^{(3)}$}};" << std::endl;
		// myfile << "\\node[main node] (0)  {0};" << std::endl;
		for (size_t i = 1; i < nBoxesPerLevel[nLevels]; i++) {
			if(i < 4)
				myfile << "\\node[main node] (";
			else
				myfile << "\\node[main node,color=white] (";
			myfile << i;
			myfile << ")   [right of=";
			myfile << i-1;
			myfile << "]  {\\scalebox{2.5}{$\\displaystyle x_{";
			myfile << i;
			myfile << "}^{(3)}$}};" << std::endl;
		}
		//self
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < 4; i++) {
			myfile << "(";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=red, loop below] node [color=red, pos=0.5, above] {} (";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;

		//P2P
		for (size_t i = 0; i < 4; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for (size_t n = 0; n < 8; n++) {//neighbors of k
				int nn = tree[nLevels][i].neighborNumbers[n];//P2P(j,nn)
				if (nn != -1) {
					if (i < nn) {//avoiding two sided
						myfile << "(";
						myfile << i;
						if (i%2 == 0) {
							myfile << ") edge[line width=2pt, draw=red, bend right, <->] node [color=red, pos=0.5, above] {} (";
						}
						else {
							myfile << ") edge[line width=2pt, draw=red, bend left, <->] node [color=red, pos=0.5, above] {} (";
						}
						myfile << nn;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		//z3
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			if(i<4)
				myfile << "\\node[main node] (z";
			else
				myfile << "\\node[main node,color=white] (z";
			myfile << i;
			myfile << ") [above of=";
			myfile << i;
			myfile << ",yshift=2cm] {\\scalebox{2.5}{$\\displaystyle z_{";
			myfile << i;
			myfile << "}^{(3)}$}};" << std::endl;
		}
		//x3<->z3
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < 4; i++) {
			myfile << "(";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=green, <->] node [color=green, pos=0.5, right] {} (z";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;
		//y3
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			if(i<4)
				myfile << "\\node[main node] (y";
			else
				myfile << "\\node[main node,color=white] (y";
			myfile << i;
			myfile << ") [above of=z";
			myfile << i;
			// myfile << ",yshift=6cm] {};";
			myfile << ",yshift=2.5cm] {\\scalebox{2.5}{$\\displaystyle y_{";
			myfile << i;
			myfile << "}^{(3)}$}};" << std::endl;
		}
		//z3<->y3
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < 4; i++) {
			myfile << "(z";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=blue, <->] node [color=green, pos=0.5, right] {} (y";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;
		//M2L3
		for (size_t i = 0; i < 4; i++) {
		// for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<16; ++in) {
				if(tree[nLevels][i].innerNumbers[in] != -1) {
					int kIL = tree[nLevels][i].innerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(y";
						myfile << i;
						if (kIL%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		for (size_t i = 0; i < 4; i++) {
		// for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<24; ++in) {
				if(tree[nLevels][i].outerNumbers[in] != -1) {
					int kIL = tree[nLevels][i].outerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(y";
						myfile << i;
						if (kIL%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		//z2
		int u=0;
		for (size_t i = 1; i < nBoxesPerLevel[nLevels]; i+=4) {
			// \node[main node, above=of $(8.north)!0.5!(9.north)$] (10) {10};
			if (i<4)
				myfile << "\\node[main node, above=of $(y";
			else
				myfile << "\\node[main node, color=white, above=of $(y";
			myfile << i;
			myfile << ".north)!0.5!(y";
			myfile << i+1;
			myfile << ".north)$,yshift=2cm] (2z";
			myfile << u;
			myfile << ") {\\scalebox{2.5}{$\\displaystyle z_{";
			myfile << u;
			myfile << "}^{(2)}$}};" << std::endl;
			u++;
		}
		// //y2<->z2
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < 4; i++) {
			myfile << "(y";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=green, <->] node [color=green, pos=0.5, right] {} (2z";
			myfile << int(i/4);
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;
		// //z2<->y2
		for (size_t i = 0; i < nBoxesPerLevel[nLevels-1]; i++) {
			if(i<1)
				myfile << "\\node[main node] (2y";
			else
				myfile << "\\node[main node, color=white] (2y";
			myfile << i;
			myfile << ") [above of=2z";
			myfile << i;
			// myfile << ",yshift=6cm] {};";
			myfile << ",yshift=2cm] {\\scalebox{2.5}{$\\displaystyle y_{";
			myfile << i;
			myfile << "}^{(2)}$}};" << std::endl;
		}
		// //z3<->y3
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < 1; i++) {
			myfile << "(2z";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=blue, <->] node [color=green, pos=0.5, right] {} (2y";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;

		//M2L2
		for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<16; ++in) {
				if(tree[nLevels-1][i].innerNumbers[in] != -1) {
					int kIL = tree[nLevels-1][i].innerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(2y";
						myfile << i;
						if (kIL%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (2y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (2y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<24; ++in) {
				if(tree[nLevels-1][i].outerNumbers[in] != -1) {
					int kIL = tree[nLevels-1][i].outerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(2y";
						myfile << i;
						if (kIL%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (2y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (2y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		myfile << "\\end{tikzpicture}" << std::endl;
		myfile << "\\end{document}" << std::endl;
		}

	void outputGraphExtendedSparseMatrix(std::string filename) {
		std::ofstream myfile;
		myfile.open(filename.c_str());
		myfile << "\\documentclass{article}" << std::endl;
		myfile << "\\usepackage{geometry}" << std::endl;
 		myfile << "\\geometry{left=20mm, top=20mm,papersize={120cm,60cm},}" << std::endl;
		myfile << "\\usepackage{pgfplots}" << std::endl;
		myfile << "\\usepackage{amsmath}" << std::endl;
		myfile << "\\usepackage{amssymb}" << std::endl;
		myfile << "\\usetikzlibrary{arrows}" << std::endl;
		myfile << "\\usetikzlibrary{shapes}" << std::endl;
		myfile << "\\usetikzlibrary{arrows.meta,calc,positioning,quotes}" << std::endl;
		myfile << "\\pgfplotsset{compat=newest,width=240cm, height=50cm,xtick pos=left, ytick pos=left,}" << std::endl;
		myfile << "\\begin{document}" << std::endl;
		myfile << "\\begin{tikzpicture}[->,>=stealth',auto,node distance=1cm, thick,main node/.style={circle,draw,minimum size=.1cm,font=\\sffamily\\Large\\bfseries}]" << std::endl;
		myfile << "\\node[main node] (0)  {};" << std::endl;
		// myfile << "\\node[main node] (0)  {0};" << std::endl;
		for (size_t i = 1; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "\\node[main node] (";
			myfile << i;
			myfile << ")   [right of=";
			myfile << i-1;
			myfile << "]  {";
			// myfile << i;
			myfile << "};" << std::endl;
		}
		//self
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "(";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=red, loop below] node [color=red, pos=0.5, above] {} (";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;

		//P2P
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for (size_t n = 0; n < 8; n++) {//neighbors of k
				int nn = tree[nLevels][i].neighborNumbers[n];//P2P(j,nn)
				if (nn != -1) {
					if (i < nn) {//avoiding two sided
						myfile << "(";
						myfile << i;
						if (i%2 == 0) {
							myfile << ") edge[line width=2pt, draw=red, bend right, <->] node [color=red, pos=0.5, above] {} (";
						}
						else {
							myfile << ") edge[line width=2pt, draw=red, bend left, <->] node [color=red, pos=0.5, above] {} (";
						}
						myfile << nn;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		//z3
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "\\node[main node] (z";
			myfile << i;
			myfile << ") [above of=";
			myfile << i;
			myfile << ",yshift=6cm] {};" << std::endl;
		}
		//x3<->z3
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "(";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=green, <->] node [color=green, pos=0.5, right] {} (z";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;
		//y3
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "\\node[main node] (y";
			myfile << i;
			myfile << ") [above of=z";
			myfile << i;
			// myfile << ",yshift=6cm] {};";
			myfile << ",yshift=6cm] {};" << std::endl;
		}
		//z3<->y3
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "(z";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=blue, <->] node [color=green, pos=0.5, right] {} (y";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;
		//M2L3
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
		// for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<16; ++in) {
				if(tree[nLevels][i].innerNumbers[in] != -1) {
					int kIL = tree[nLevels][i].innerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(y";
						myfile << i;
						if (i%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
		// for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<24; ++in) {
				if(tree[nLevels][i].outerNumbers[in] != -1) {
					int kIL = tree[nLevels][i].outerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(y";
						myfile << i;
						if (i%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		//z2
		int u=0;
		for (size_t i = 1; i < nBoxesPerLevel[nLevels]; i+=4) {
			// \node[main node, above=of $(8.north)!0.5!(9.north)$] (10) {10};
			myfile << "\\node[main node, above=of $(y";
			myfile << i;
			myfile << ".north)!0.5!(y";
			myfile << i+1;
			myfile << ".north)$,yshift=6cm] (2z";
			myfile << u;
			myfile << ") {};" << std::endl;
			u++;
		}
		//y2<->z2
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < nBoxesPerLevel[nLevels]; i++) {
			myfile << "(y";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=green, <->] node [color=green, pos=0.5, right] {} (2z";
			myfile << int(i/4);
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;
		//z2<->y2
		for (size_t i = 0; i < nBoxesPerLevel[nLevels-1]; i++) {
			myfile << "\\node[main node] (2y";
			myfile << i;
			myfile << ") [above of=2z";
			myfile << i;
			// myfile << ",yshift=6cm] {};";
			myfile << ",yshift=6cm] {};" << std::endl;
		}
		//z3<->y3
		myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
		for (size_t i = 0; i < nBoxesPerLevel[nLevels-1]; i++) {
			myfile << "(2z";
			myfile << i;
			myfile << ") edge[line width=2pt, draw=blue, <->] node [color=green, pos=0.5, right] {} (2y";
			myfile << i;
			myfile << ")" << std::endl;
		}
		myfile << ";" << std::endl;

		//M2L2
		for (size_t i = 0; i < nBoxesPerLevel[nLevels-1]; i++) {
		// for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<16; ++in) {
				if(tree[nLevels-1][i].innerNumbers[in] != -1) {
					int kIL = tree[nLevels-1][i].innerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(2y";
						myfile << i;
						if (i%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (2y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (2y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		for (size_t i = 0; i < nBoxesPerLevel[nLevels-1]; i++) {
		// for (size_t i = 0; i < 1; i++) {
			myfile << "\\path[every node/.style={font=\\sffamily\\small}]" << std::endl;
			for(int in=0; in<24; ++in) {
				if(tree[nLevels-1][i].outerNumbers[in] != -1) {
					int kIL = tree[nLevels-1][i].outerNumbers[in];
					if (i < kIL) {//avoiding two sided
						myfile << "(2y";
						myfile << i;
						if (i%2 == 0) {
							myfile << ") edge[line width=2pt, draw=orange, bend right, <->] node [color=red, pos=0.5, above] {} (2y";
						}
						else {
							myfile << ") edge[line width=2pt, draw=orange, bend left, <->] node [color=red, pos=0.5, above] {} (2y";
						}
						myfile << kIL;
						myfile << ")" << std::endl;
					}
				}
			}
			myfile << ";" << std::endl;
		}
		myfile << "\\end{tikzpicture}" << std::endl;
		myfile << "\\end{document}" << std::endl;
		}

	void outputExtendedSparseMatrix(std::string filename) {
		double xcenter = 0.0;
		double ycenter = 0.0;
		double Lx = L;
		double Ly = L;

		std::ofstream myfile;
		myfile.open(filename.c_str());
		myfile << "\\documentclass{standalone}" << std::endl;
		myfile << "\\usepackage{tikz}" << std::endl;
		myfile << "\\begin{document}" << std::endl;
		myfile << "\\begin{tikzpicture}" << std::endl;
		// x^2
		int j=nLevels;
		int sizeOfMatrix = 3*nBoxesPerLevel[nLevels];
		for (size_t i = nLevels-1; i >= 2; i--) {
			sizeOfMatrix += 2*nBoxesPerLevel[i];
		}
		std::cout << "sizeOfMatrix: " << sizeOfMatrix << std::endl;
		double x,x2;
		double y;
		double diameter = 10.0/sizeOfMatrix;
		// myfile << "\\draw[step=";
		// myfile << diameter;
		// myfile << "] (0,0) grid +(10,10);" << std::endl;
		myfile << "\\draw (0,0) rectangle +(10,10);" << std::endl;

		// x3
		for (int i=0; i<nBoxesPerLevel[nLevels]; ++i) {//rows
			//self
			x = 0.0;
			y = 10.0-diameter;
			for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
				myfile << "\\draw[draw=none,fill=red] (" << x << ",";
				myfile << y << ") rectangle +(";
				myfile << diameter << ",";
				myfile << diameter << ");" << std::endl;
				x += 2*diameter;
				y -= 2*diameter;
			}
			//neighbors
			for (size_t n = 0; n < 8; n++) {//neighbors of k
				int nn = tree[j][i].neighborNumbers[n];//P2P(j,nn)
				if (nn != -1) {
					double diameter = 10.0/sizeOfMatrix;
					double kx = i*2*diameter;
					double ky = 10-(i+1)*2*diameter;
					double nnx = nn*2*diameter;
					double nny = 10-(nn+1)*2*diameter+diameter;
					myfile << "\\draw[draw=none,fill=red] (" << kx << ",";
					myfile << nny << ") rectangle +(";
					myfile << diameter << ",";
					myfile << diameter << ");" << std::endl;
				}
			}
		}
			//L2P
			x = diameter;
			y = 10.0-diameter;
			for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
				myfile << "\\draw[draw=none,fill=green] (" << x << ",";
				myfile << y << ") rectangle +(";
				myfile << diameter << ",";
				myfile << diameter << ");" << std::endl;
				x += 2*diameter;
				y -= 2*diameter;
			}
			//z3
			// for (int i=0; i<nBoxesPerLevel[nLevels]; ++i) {//rows
				//V
				x = 0.0;
				x2 = 2*nBoxesPerLevel[nLevels]*diameter;
				y = 10.0-2*diameter;
				for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
					myfile << "\\draw[draw=none,fill=green] (" << x << ",";
					myfile << y << ") rectangle +(";
					myfile << diameter << ",";
					myfile << diameter << ");" << std::endl;
					x += 2*diameter;
					//-I
					myfile << "\\draw[draw=none,fill=blue] (" << x2 << ",";
					myfile << y << ") rectangle +(";
					myfile << diameter << ",";
					myfile << diameter << ");" << std::endl;
					if ((k+1)%4 == 0)
						x2 += 2*diameter;
					else
					x2 += diameter;
					y -= 2*diameter;
				}
			//y3 = x2
			x = diameter;
			y = 10.0-2*diameter*nBoxesPerLevel[nLevels]-diameter;
			for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
				//-I
				myfile << "\\draw[draw=none,fill=blue] (" << x << ",";
				myfile << y << ") rectangle +(";
				myfile << diameter << ",";
				myfile << diameter << ");" << std::endl;
				x += 2*diameter;
				if ((k+1)%4 == 0)
					y -= 2*diameter;
				else
				y -= diameter;
			}
			for (int k=0; k<nBoxesPerLevel[nLevels]; ++k) {
				//IL
				for(int in=0; in<16; ++in) {
					if(tree[j][k].innerNumbers[in] != -1) {
						int kIL = tree[j][k].innerNumbers[in];
						double kx = kIL*diameter + 2*diameter*nBoxesPerLevel[nLevels] + int(kIL/4)*diameter;
						double nny = 10-((k+1)*diameter + 2*diameter*nBoxesPerLevel[nLevels]) - int(k/4)*diameter;
						myfile << "\\draw[draw=none,fill=orange] (" << kx << ",";
						myfile << nny << ") rectangle +(";
						myfile << diameter << ",";
						myfile << diameter << ");" << std::endl;
					}
				}
				for(int on=0; on<24; ++on) {
					if(tree[j][k].outerNumbers[on] != -1) {
						int kIL = tree[j][k].outerNumbers[on];
						double kx = kIL*diameter + 2*diameter*nBoxesPerLevel[nLevels] + int(kIL/4)*diameter;
						double nny = 10-((k+1)*diameter+ 2*diameter*nBoxesPerLevel[nLevels]) - int(k/4)*diameter;
						myfile << "\\draw[draw=none,fill=orange] (" << kx << ",";
						myfile << nny << ") rectangle +(";
						myfile << diameter << ",";
						myfile << diameter << ");" << std::endl;
					}
				}
			}
			// z2
			x = 2*diameter*nBoxesPerLevel[nLevels];
			y = 10-(4*diameter+ 2*diameter*nBoxesPerLevel[nLevels])-diameter;
			x2 = 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1];
			double y2 = 10-(4*diameter+ 2*diameter*nBoxesPerLevel[nLevels])-diameter;
			for (int k=0; k<nBoxesPerLevel[nLevels-1]; ++k) {
				myfile << "\\draw[draw=none,fill=green] (" << x << ",";
				myfile << y << ") rectangle +(";
				myfile << 4*diameter << ",";
				myfile << diameter << ");" << std::endl;
				x += 5*diameter;
				y -= 5*diameter;
				//-I
				myfile << "\\draw[draw=none,fill=blue] (" << x2 << ",";
				myfile << y2 << ") rectangle +(";
				myfile << diameter << ",";
				myfile << diameter << ");" << std::endl;
				x2 += diameter;
				y2 -= 5*diameter;
			}
			y = 10-(2*diameter*nBoxesPerLevel[nLevels])-4*diameter;
			x = 2*diameter*nBoxesPerLevel[nLevels]+4*diameter;
			for (int k=0; k<nBoxesPerLevel[nLevels-1]; ++k) {
				myfile << "\\draw[draw=none,fill=green] (" << x << ",";
				myfile << y << ") rectangle +(";
				myfile << diameter << ",";
				myfile << 4*diameter << ");" << std::endl;
				x += 5*diameter;
				y -= 5*diameter;
			}
			//y2
			j = nLevels-1;
			x = 2*diameter*nBoxesPerLevel[nLevels] + 4*diameter;
			for (int k=0; k<nBoxesPerLevel[nLevels-1]; ++k) {
				double nny = 10-((k+1)*diameter + 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1]);
				myfile << "\\draw[draw=none,fill=blue] (" << x << ",";
				myfile << nny << ") rectangle +(";
				myfile << diameter << ",";
				myfile << diameter << ");" << std::endl;
				x += 4*diameter;
			}
			for (int k=0; k<nBoxesPerLevel[nLevels-1]; ++k) {
				//IL
				for(int in=0; in<16; ++in) {
					if(tree[j][k].innerNumbers[in] != -1) {
						int kIL = tree[j][k].innerNumbers[in];
						double kx = kIL*diameter + 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1];
						double nny = 10-((k+1)*diameter + 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1]);
						myfile << "\\draw[draw=none,fill=orange] (" << kx << ",";
						myfile << nny << ") rectangle +(";
						myfile << diameter << ",";
						myfile << diameter << ");" << std::endl;
					}
				}
				for(int on=0; on<24; ++on) {
					if(tree[j][k].outerNumbers[on] != -1) {
						int kIL = tree[j][k].outerNumbers[on];
						double kx = kIL*diameter + 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1];
						double nny = 10-((k+1)*diameter + 3*diameter*nBoxesPerLevel[nLevels] + diameter*nBoxesPerLevel[nLevels-1]);
						myfile << "\\draw[draw=none,fill=orange] (" << kx << ",";
						myfile << nny << ") rectangle +(";
						myfile << diameter << ",";
						myfile << diameter << ");" << std::endl;
					}
				}
			}

		myfile << "\\end{tikzpicture}" << std::endl;
		myfile << "\\end{document}" << std::endl;
		myfile.close();
	}

};

#endif
