#ifndef _AFMM2DTree_HPP__
#define _AFMM2DTree_HPP__
class AFMM2DBox {
public:
	int boxNumber;
	int parentNumber;
	int childrenNumbers[4];
	int neighborNumbers[8];
	int innerNumbers[16];
	int outerNumbers[24];

	AFMM2DBox () {
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
	}
	pts2D center;
	Vec charges;

	std::map<int, Mat> M2L;
	Mat L2P, P2M;					//	Transfer from multipoles of 4 children to multipoles of parent.
	//	The following will be stored only at the leaf nodes
	std::vector<pts2D> chebNodes;
	std::vector<pts2D> particle_loc;
  std::vector<int> chargeLocations;

	Vec outgoing_charges;//equivalent densities {f_{k}^{B,o}}
	Vec incoming_charges;//equivalent densities {f_{k}^{B,i}}
	Vec incoming_potential;//check potentials {u_{k}^{B,i}}
	std::vector<int> incoming_chargePoints;//equivalent points {y_{k}^{B,i}}
	std::vector<int> incoming_checkPoints;//check points {x_{k}^{B,i}}
	std::vector<int> outgoing_checkPoints;//equivalent points {y_{k}^{B,i}}
	std::vector<int> outgoing_chargePoints;//check points {x_{k}^{B,i}}
	Vec potential;//used only at leaf level
};

template <typename kerneltype>
class AFMM2DTree {
public:
	kerneltype* K;
	int nLevels;			//	Number of levels in the tree.
	int N;					//	Number of particles.
	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<std::vector<AFMM2DBox> > tree;	//	The tree storing all the information.
	std::vector<double> boxRadius;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	double L;				//	Semi-length of the simulation box.
	double ACA_epsilon;
  std::vector<pts2D> gridPoints; //all particles in domain
	int TOL_POW;

// public:
	AFMM2DTree(kerneltype* K, double L, int N, int nParticlesInLeafAlong1D, int TOL_POW, std::vector<pts2D>& gridPoints) {
		this->K					=	K;
		this->L					=	L;
		this->nLevels		=	ceil(log(double(N)/nParticlesInLeafAlong1D/nParticlesInLeafAlong1D)/log(4));
		this->TOL_POW = TOL_POW;
		this->gridPoints = gridPoints;
		nBoxesPerLevel.push_back(1);
		boxRadius.push_back(L);
		for (int k=1; k<=nLevels; ++k) {
			nBoxesPerLevel.push_back(4*nBoxesPerLevel[k-1]);
			boxRadius.push_back(0.5*boxRadius[k-1]);
		}
		this->N					=	N;//nParticlesInLeaf*nBoxesPerLevel[nLevels];
		this->ACA_epsilon = 1.0e-10;
	}

	void createTree() {
		//	First create root and add to tree
		AFMM2DBox root;
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
		std::vector<AFMM2DBox> rootLevel;
		rootLevel.push_back(root);
		tree.push_back(rootLevel);

		for (int j=1; j<=nLevels; ++j) {
			std::vector<AFMM2DBox> level;
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				AFMM2DBox box;
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

	void assignChargeLocations() {
		// K->particles = Nodes;//object of base class FMM_Matrix
		for (size_t i = 0; i < N; i++) {
			tree[0][0].chargeLocations.push_back(i);
		}
		for (size_t j = 0; j < nLevels; j++) { //assign particles to its children
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				int J = j+1;
				int Kp = 4*k;
				for (size_t i = 0; i < tree[j][k].chargeLocations.size(); i++) {
					int index = tree[j][k].chargeLocations[i];
					if (K->particles_X[index].x <= tree[j][k].center.x) { //children 0,3
						if (K->particles_X[index].y <= tree[j][k].center.y) { //child 0
							tree[J][Kp].chargeLocations.push_back(index);
						}
						else { //child 3
							tree[J][Kp+3].chargeLocations.push_back(index);
						}
					}
					else { //children 1,2
						if (K->particles_X[index].y <= tree[j][k].center.y) { //child 1
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

	void assignNonLeafChargeLocations() {
		for (int j = nLevels-1; j >= 1; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				tree[j][k].chargeLocations.clear();
				for (size_t c = 0; c < 4; c++) {
					tree[j][k].chargeLocations.insert(tree[j][k].chargeLocations.end(), tree[j+1][4*k+c].chargeLocations.begin(), tree[j+1][4*k+c].chargeLocations.end());
				}
				// std::cout << "tree[j][k].charges: " << tree[j][k].charges.size() << std::endl;
			}
		}
	}


	// void assignLeafCharges(Vec &chargesVec) {
	// 	int index = 0;
	// 	for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
	// 		int boxNum = boxNumbers[k];
	// 		int NumParticles = NumberOfParticlesInLeaves[k];
	// 		tree[nLevels][boxNum].charges = Vec(NumParticles);
	// 		for (size_t i = 0; i < NumParticles; i++) {
	// 			tree[nLevels][boxNum].charges(i) = chargesVec(index+i);
	// 		}
	// 		index += NumParticles;
	// 	}
	// }

	void assignLeafCharges(Vec &charges) {
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			tree[nLevels][k].charges	=	Eigen::VectorXd::Zero(tree[nLevels][k].chargeLocations.size());
			for (size_t i = 0; i < tree[nLevels][k].chargeLocations.size(); i++) {
				int index = tree[nLevels][k].chargeLocations[i];
				tree[nLevels][k].charges[i]	=	charges[index];
			}
		}
	}


	void getNodes() {
		for (int j=nLevels; j>=2; j--) {
			getNodes_incoming_level(j);
			getNodes_outgoing_level(j);
		}
	}

	void getNodes_incoming_level(int j) {
		int ComputedRank;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			// std::cout << "Ij: " << j << "	k: " << k << std::endl;
			getNodes_incoming_box(j, k, ComputedRank);
		}
	}

	void getNodes_outgoing_level(int j) {
		int ComputedRank;
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			// std::cout << "Oj: " << j << "	k: " << k << std::endl;
			getNodes_outgoing_box(j, k, ComputedRank);
		}
	}

	void getParticlesFromChildren_row(int j, int k, std::vector<int>& searchNodes) {
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

	void getParticlesFromChildren_col(int j, int k, std::vector<int>& searchNodes) {
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

	void getNodes_incoming_box(int j, int k, int& ComputedRank) {
		std::vector<int> boxA_Nodes;
		getParticlesFromChildren_row(j, k, boxA_Nodes);
		std::vector<int> IL_Nodes;//indices
		for(int in=0; in<16; ++in) {
			if(tree[j][k].innerNumbers[in] != -1) {
				int kIL = tree[j][k].innerNumbers[in];
				std::vector<int> chargeLocations;
				getParticlesFromChildren_col(j, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}
		}
		for(int on=0; on<24; ++on) {
			if(tree[j][k].outerNumbers[on] != -1) {
				int kIL = tree[j][k].outerNumbers[on];
				std::vector<int> chargeLocations;
				getParticlesFromChildren_col(j, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}
		}
		if(IL_Nodes.size() > 0 && boxA_Nodes.size() > 0) {
			std::vector<int> row_bases, col_bases;
			Mat Ac, Ar, L, R;
			LowRank* LR		=	new LowRank(K, pow(10,-TOL_POW), boxA_Nodes, IL_Nodes);
			// std::cout << "before " << std::endl;
			LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, Ar, L, R);
			// std::cout << "after " << std::endl;
			// std::cout << "k: " << k << "	CR: " << ComputedRank << std::endl;
			// LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, Ar);
			if(ComputedRank > 0) {
				for (int r = 0; r < row_bases.size(); r++) {
					tree[j][k].incoming_checkPoints.push_back(boxA_Nodes[row_bases[r]]);
				}
				for (int c = 0; c < col_bases.size(); c++) {
					tree[j][k].incoming_chargePoints.push_back(IL_Nodes[col_bases[c]]);
				}

				Mat temp = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheRight>(Ac);
				tree[j][k].L2P = L.triangularView<Eigen::Lower>().solve<Eigen::OnTheRight>(temp);

				// Mat D = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][k].incoming_chargePoints);
				// Eigen::PartialPivLU<Mat> D_T = Eigen::PartialPivLU<Mat>(D.transpose());
				// Eigen::PartialPivLU<Mat> D_T_T = Eigen::PartialPivLU<Mat>(D);
				// Mat Error_LR_D = L*R - D;
				// std::cout << "Ik-----------------: " << k << std::endl;
				// std::cout << "E_LR_D: " << Error_LR_D.norm()/D.norm() << std::endl;
				// // std::cout << "mat: " << std::endl << K->getMatrix(IL_Nodes,boxA_Nodes) << std::endl;
				// Mat Error_AcAr_LR = K->getMatrix(boxA_Nodes,IL_Nodes) - Ac*D_T_T.solve(Ar);//tree[j][k].L2P*Ar;//K->getMatrix(boxA_Nodes, IL_Nodes);
				// std::cout << "E: " << Error_AcAr_LR.norm() << std::endl;
				// std::cout << "ComputedRank: " << ComputedRank << std::endl;

			}
			else {
				tree[j][k].L2P = Mat::Zero(Ac.rows(),ComputedRank);
			}
		}
		else {
			tree[j][k].L2P = Mat::Zero(boxA_Nodes.size(),0);
		}
	}

	void getNodes_outgoing_box(int j, int k, int& ComputedRank) {
		std::vector<int> boxA_Nodes;
		getParticlesFromChildren_col(j, k, boxA_Nodes);
		std::vector<int> IL_Nodes;//indices
		for(int in=0; in<16; ++in) {
			if(tree[j][k].innerNumbers[in] != -1) {
				int kIL = tree[j][k].innerNumbers[in];
				std::vector<int> chargeLocations;
				getParticlesFromChildren_row(j, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}
		}
		for(int on=0; on<24; ++on) {
			if(tree[j][k].outerNumbers[on] != -1) {
				int kIL = tree[j][k].outerNumbers[on];
				std::vector<int> chargeLocations;
				getParticlesFromChildren_row(j, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}
		}
		if(IL_Nodes.size() > 0 && boxA_Nodes.size() > 0) {
			std::vector<int> row_bases, col_bases;
			Mat Ac, Ar, L, R;
			LowRank* LR		=	new LowRank(K, pow(10,-TOL_POW), IL_Nodes, boxA_Nodes);
			// std::cout << "size: " << IL_Nodes.size() << "," << boxA_Nodes.size() << std::endl;
			LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, Ar, R, L);
			// std::cout << "after " << std::endl;
			// std::cout << "k: " << k << "	CR: " << ComputedRank << std::endl;
			// LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, Ar);
			if(ComputedRank > 0) {
				for (int r = 0; r < row_bases.size(); r++) {
					tree[j][k].outgoing_checkPoints.push_back(IL_Nodes[row_bases[r]]);
				}
				for (int c = 0; c < col_bases.size(); c++) {
					tree[j][k].outgoing_chargePoints.push_back(boxA_Nodes[col_bases[c]]);
				}

				Mat temp = L.triangularView<Eigen::Lower>().solve<Eigen::OnTheLeft>(Ar);
				tree[j][k].P2M = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(temp);

				// if(k==0 || k==5){
				// Mat D = K->getMatrix(tree[j][k].outgoing_checkPoints, tree[j][k].outgoing_chargePoints);
				// Eigen::PartialPivLU<Mat> D_T = Eigen::PartialPivLU<Mat>(D.transpose());
				// Eigen::PartialPivLU<Mat> D_T_T = Eigen::PartialPivLU<Mat>(D);
				// Mat Error_LR_D = L*R - D;
				// std::cout << "Ok-----------------: " << k << std::endl;
				// std::cout << "E_LR_D: " << Error_LR_D.norm()/D.norm() << std::endl;
				// std::cout << "mat: " << std::endl << K->getMatrix(IL_Nodes,boxA_Nodes) << std::endl;
				// std::cout << "Ac*D_T_T.solve(Ar): " << std::endl << Ac*D_T_T.solve(Ar) << std::endl;
				// std::cout << "Err: " << std::endl << K->getMatrix(IL_Nodes,boxA_Nodes) - Ac*D_T_T.solve(Ar) << std::endl;
			// 	Mat Error_AcAr_LR = K->getMatrix(IL_Nodes,boxA_Nodes) - Ac*D_T_T.solve(Ar);//tree[j][k].L2P*Ar;//K->getMatrix(boxA_Nodes, IL_Nodes);
			// 	std::cout << "E: " << Error_AcAr_LR.norm() << std::endl;
			// 	std::cout << "ComputedRank: " << ComputedRank << std::endl;
			// // }
			// 	std::cout << "L: " << std::endl << L << std::endl;
			// 	std::cout << "R: " << std::endl << R << std::endl;
			// 	std::cout << "LR: " << std::endl << L*R << std::endl;
			// 	std::cout << "RL: " << std::endl << R*L << std::endl;
			// 	std::cout << "D: " << std::endl << D << std::endl;
			// 	std::cout << "LR-D: " << std::endl << L*R-D << std::endl;
			// 	exit(0);

			}
			else {
				tree[j][k].P2M = Mat::Zero(ComputedRank,Ar.cols());
			}
		}
		else {
			tree[j][k].P2M = Mat::Zero(0,boxA_Nodes.size());
		}
	}


	void assemble_M2L() {
		#pragma omp parallel for
		for (size_t j = 2; j <= nLevels; j++) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				#pragma omp parallel for
				for(int in=0; in<16; ++in) {
					if(tree[j][k].innerNumbers[in] != -1) {
						int kIL = tree[j][k].innerNumbers[in];
						// if (tree[j][k].M2L[kIL].size() == 0)
							// tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].incoming_checkPoints);
							tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].outgoing_chargePoints);
						// if (tree[j][kIL].M2L[k].size() == 0)
						// 	tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose();
					}
				}
				#pragma omp parallel for
				for(int on=0; on<24; ++on) {
					if(tree[j][k].outerNumbers[on] != -1) {
						int kIL = tree[j][k].outerNumbers[on];
						// if (tree[j][k].M2L[kIL].size() == 0)
							// tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].incoming_checkPoints);
							tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].outgoing_chargePoints);
						// if (tree[j][kIL].M2L[k].size() == 0)
						// 	tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose();
					}
				}
			}
		}
	}

	void evaluate_M2M() {
		// std::cout << "nLevels: " << nLevels << std::endl;
		for (size_t j = nLevels; j >= 2; j--) {
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				Vec source_densities;
				if (j==nLevels) {
					source_densities = tree[j][k].charges;
				}
				else {
					int J = j+1;
					int Veclength = 0;
					for (int child = 0; child < 4; child++) {
						Veclength += tree[J][4*k+child].outgoing_charges.size();
					}
					source_densities = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
					int start = 0;
					for (int child = 0; child < 4; child++) {
						int NumElem = tree[J][4*k+child].outgoing_charges.size();
						source_densities.segment(start, NumElem) = tree[J][4*k+child].outgoing_charges;
						start += NumElem;
					}
				}
				// std::cout << "k: " << k << std::endl;
				// std::cout << "P2M: " << tree[j][k].P2M.rows() << "," << tree[j][k].P2M.cols() << std::endl;
				// std::cout << "sd: " << source_densities.size() << std::endl;
				// tree[j][k].outgoing_charges = tree[j][k].L2P.transpose()*source_densities;//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
				tree[j][k].outgoing_charges = tree[j][k].P2M*source_densities;//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
			}
		}
	}

	void evaluate_M2L() {
		#pragma omp parallel for
		for (int j=2; j<=nLevels; ++j) {
			// #pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {//BoxA
				tree[j][k].incoming_potential	=	Vec::Zero(tree[j][k].incoming_checkPoints.size());
				// #pragma omp parallel for
				for(int in=0; in<16; ++in) {
					if(tree[j][k].innerNumbers[in] != -1) {
						int kIL = tree[j][k].innerNumbers[in];
						tree[j][k].incoming_potential += tree[j][k].M2L[kIL]*tree[j][kIL].outgoing_charges;
					}
				}
				// #pragma omp parallel for
				for(int on=0; on<24; ++on) {
					if(tree[j][k].outerNumbers[on] != -1) {
						int kIL = tree[j][k].outerNumbers[on];
						tree[j][k].incoming_potential += tree[j][k].M2L[kIL]*tree[j][kIL].outgoing_charges;
					}
				}
			}
		}
	}

	void evaluate_L2L() {
		for (size_t j = 2; j <= nLevels; j++) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				if (j != nLevels) {
					Vec temp = tree[j][k].L2P*tree[j][k].incoming_potential;
					int start = 0;
					for (size_t c = 0; c < 4; c++) {
						tree[j+1][4*k+c].incoming_potential += temp.segment(start, tree[j+1][4*k+c].incoming_checkPoints.size());
						start += tree[j+1][4*k+c].incoming_checkPoints.size();
					}
				}
				else {
					tree[j][k].potential = tree[j][k].L2P*tree[j][k].incoming_potential;
				}
			}
		}
	}

	void evaluate_NearField() { // evaluating at chargeLocations
		#pragma omp parallel for
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			#pragma omp parallel for
			for (size_t n = 0; n < 8; n++) {
				int nn = tree[nLevels][k].neighborNumbers[n];
				if(nn != -1) {
					Mat R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][nn].chargeLocations);
					tree[nLevels][k].potential += R*tree[nLevels][nn].charges;
				}
			}
			Mat R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][k].chargeLocations);
			tree[nLevels][k].potential += R*tree[nLevels][k].charges; //self Interaction
		}
	}

	void collectPotential(Vec &potential) {
		potential = Vec::Zero(N);
		int start = 0;
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			potential.segment(start, tree[nLevels][k].potential.size()) = tree[nLevels][k].potential;
			start += tree[nLevels][k].potential.size();
		}
	}

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

	int getMaxRank() {
		int sum = 0;
		int NumBoxes = 0;
		// #pragma omp parallel for
		for (size_t j = 2; j <= nLevels; j++) {
			// NumBoxes += nBoxesPerLevel[j];
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				// NumBoxes += 1;
				if(sum < tree[j][k].incoming_checkPoints.size()){
					sum = tree[j][k].incoming_checkPoints.size();
				}
			}
		}
		return sum;
		// return sum/NumBoxes;
	}

	double compute_error(int nBox) { // evaluating at chargeLocations
		Vec true_potential = Vec::Zero(tree[nLevels][nBox].chargeLocations.size());
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			Mat R = K->getMatrix(tree[nLevels][nBox].chargeLocations, tree[nLevels][k].chargeLocations);
			true_potential += R*tree[nLevels][k].charges;
		}
		Vec errVec = true_potential - tree[nLevels][nBox].potential;
		double error = errVec.norm()/true_potential.norm();
		return error;
	}
};

#endif
