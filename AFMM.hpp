#include "AFMM2DTree.hpp"

template <typename kerneltype>
class AFMM {
public:
	AFMM2DTree<kerneltype>* A;
		AFMM(int nParticlesInLeafAlong1D, double L, int TOL_POW, std::vector<pts2D> gridPoints, kerneltype* mykernel) {
			unsigned n_Dimension    =       2;  //      Dimension.

			unsigned count_Location =       0;
			int N = gridPoints.size();

			A	=	new AFMM2DTree<kerneltype>(mykernel, L, N, nParticlesInLeafAlong1D, TOL_POW, gridPoints);

			A->createTree();
			A->assign_Tree_Interactions();
			A->assign_Center_Location();
			// A->assignLeafChargeLocations();

			A->assignChargeLocations();
			A->assignNonLeafChargeLocations();

			A->getNodes();
			A->assemble_M2L();
		}

		void MatVecProduct(Vec &chargesVec, Vec &AFMM_Ab) {
			A->assignLeafCharges(chargesVec);
			A->evaluate_M2M();
			A->evaluate_M2L();
			A->evaluate_L2L();
			A->evaluate_NearField();
			A->collectPotential(AFMM_Ab);
			A->reorder(AFMM_Ab);
		}
		~AFMM() {
			delete A;
		};
};
