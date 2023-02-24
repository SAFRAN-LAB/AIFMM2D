#include "AIFMMTree.hpp"

template <typename kerneltype>
class AIFMM {
public:
  FMM2DTree<kerneltype>* A;
  AIFMM(kerneltype* mykernel, int N, int nLevels, int TOL_POW, double* locations) {
    A	=	new FMM2DTree<kerneltype>(mykernel, int(N), int(nLevels), TOL_POW, locations);
    A->createTree();
  	A->assign_Tree_Interactions();
    A->assign_Center_Location();
  	A->assignChargeLocations();
    A->assignNonLeafChargeLocations();//actually it doesnt assign; clears it
    A->getNodes();
    A->assemble_M2L();
    A->initialise_phase();
    A->initialise_P2P_Leaf_Level();
  }

  void factorize() {
  	A->eliminate_phase_efficient();
  }
  Vec solve(Vec &rhs) {
    A->assign_Leaf_rhs(rhs);
    A->rhs_eliminate_phase_efficient();
    A->back_substitution_phase();
    Vec phi;
    A->getx(phi);
    return phi;
  }
  void backSubstitute(Vec &rhs) {
    A->assign_Leaf_rhs(rhs);
    A->rhs_eliminate_phase_efficient();
    A->back_substitution_phase();
  }

  void backSubstitute1(Vec &rhs) {
    A->assign_Leaf_rhs(rhs);
  }

  void backSubstitute2() {
    A->rhs_eliminate_phase_efficient();
  }

  void backSubstitute3() {
    A->back_substitution_phase();
  }
  void getSolution(Vec &phi) {
    A->getx(phi);
    A->reorder(phi);
  }
  double getError(Vec &rhs) {
    A->assign_Leaf_rhs(rhs);
    return A->error_check();
  }
  ~AIFMM() {
    delete A;
  };
};
