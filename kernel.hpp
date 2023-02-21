//
//  kernel.hpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#ifndef __kernel_hpp__
#define __kernel_hpp__
#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

const std::complex<double> I(0.0, 1.0);
#ifdef USE_LS
  using dtype=std::complex<double>;
  using dtype_base=double;
  using Mat=Eigen::MatrixXcd;
  using Vec=Eigen::VectorXcd;
  #include "LS/LS_Assembly.hpp" //Lippmann-Schwinger kernel
#endif

#ifdef USE_Hankel
  struct pts2D {
  	double x,y;
  };
  const double PI	=	3.1415926535897932384;
  using dtype=std::complex<double>;
  using dtype_base=double;
  using Mat=Eigen::MatrixXcd;
  using Vec=Eigen::VectorXcd;
  #include <boost/math/special_functions/bessel.hpp>
  double besselJ(int n, double x) {
  	if (n >= 0) {
  		double temp = boost::math::cyl_bessel_j(double(n), x);
  		return temp;
  	}
  	else {
  		double temp = boost::math::cyl_bessel_j(double(-n), x);
  		if (-n%2 == 0)
  			return temp;
  		else
  			return -temp;
  	}
  }

  double besselY(int n, double x) {
  	if (n >= 0) {
  		double temp = boost::math::cyl_neumann(double(n), x);
  		return temp;
  	}
  	else {
  		double temp = boost::math::cyl_neumann(double(-n), x);
  		if (-n%2 == 0)
  			return temp;
  		else
  			return -temp;
  	}
  }
#endif

    // const std::complex<double> I(0.0, 1.0);
#ifdef USE_real
struct pts2D {
	double x,y;
};
const double PI	=	3.1415926535897932384;
  using dtype=double;
  using dtype_base=double;
  using Mat=Eigen::MatrixXd;
  using Vec=Eigen::VectorXd;
#endif

#ifdef USE_integralEqn
struct pts2D {
	double x,y;
};
const double PI	=	3.1415926535897932384;
  using dtype=double;
  using dtype_base=double;
  using Mat=Eigen::MatrixXd;
  using Vec=Eigen::VectorXd;
#endif

// #define EIGEN_DONT_PARALLELIZE
using namespace Eigen;

// const double PI	=	3.1415926535897932384;
#include <map>
// struct pts2D {
// 	double x,y;
// };
#include "Integral2D.hpp"

class kernel {
public:
  bool isTrans;		//	Checks if the kernel is translation invariant, i.e., the kernel is K(r).
	bool isHomog;		//	Checks if the kernel is homogeneous, i.e., K(r) = r^{alpha}.
	bool isLogHomog;	//	Checks if the kernel is log-homogeneous, i.e., K(r) = log(r^{alpha}).
	double alpha;		//	Degree of homogeneity of the kernel.
  double a;

  std::vector<pts2D> particles_X;
	std::vector<pts2D> particles_Y;

  kernel() {
	}

	virtual dtype getMatrixEntry(const unsigned i, const unsigned j) {
		std::cout << "virtual getInteraction" << std::endl;
		return 0.0;
	}

	Vec getRow(const int j, std::vector<int> col_indices) {
		int n_cols = col_indices.size();
		Vec row(n_cols);
    #pragma omp parallel for
    for(int k = 0; k < n_cols; k++) {
        row(k) = this->getMatrixEntry(j, col_indices[k]);
    }
    return row;
  }

  Vec getCol(const int k, std::vector<int> row_indices) {
		int n_rows = row_indices.size();
    Vec col(n_rows);
    #pragma omp parallel for
    for (int j=0; j<n_rows; ++j) {
			col(j) = this->getMatrixEntry(row_indices[j], k);
    }
    return col;
  }

  Vec getCol(const int n, const int k) {
    Vec col(n);
    // #pragma omp parallel for
    for (int j=0; j<n; ++j) {
			col(j) = this->getMatrixEntry(j, k);
    }
    return col;
  }

  Mat getMatrix(std::vector<int> row_indices, std::vector<int> col_indices) {
		int n_rows = row_indices.size();
		int n_cols = col_indices.size();
    Mat mat(n_rows, n_cols);
    #pragma omp parallel for
    for (int j=0; j < n_rows; ++j) {
        #pragma omp parallel for
        for (int k=0; k < n_cols; ++k) {
            mat(j,k) = this->getMatrixEntry(row_indices[j], col_indices[k]);
        }
    }
    return mat;
  }
  ~kernel() {};
};

class userkernel: public kernel {
private:
  #ifdef yestoLS
    int nCones_LFR;//number of cones to be used in DAFMM
    double TOL_POW;//currently not in use
    int nLevelsUniform; //number of levels to be constructed in case of uniform tree
    int yes2DFMM; //DFMM or FMM
  #endif
public:
  double Kii;
  double h2;
  int N;
  std::vector<pts2D> gridPoints; // location of particles in the domain
  Vec rhs;
  int nLevelsResolve;
  #ifdef yestoLS
    LS2DTree* F;
    userkernel(int nChebNodes, int treeAdaptivity, int degreeOfBases, double L, int nLevelsUniform) : kernel() {
      isTrans		=	true;
      isHomog		=	true;
      isLogHomog	=	false;
      alpha		=	-1.0;

      this->nCones_LFR			=	16;//number of cones to be used in DAFMM
      this->TOL_POW         = 9;//currently not in use
      this->nLevelsUniform  = nLevelsUniform; //number of levels to be constructed in case of uniform tree
      this->yes2DFMM				=	0; //DFMM or FMM
      this->F	=	new LS2DTree(nCones_LFR, nChebNodes, L, yes2DFMM, TOL_POW, kappa, degreeOfBases, treeAdaptivity, nLevelsUniform);
      this->nLevelsResolve = F->nLevels;

      this->gridPoints = F->gridPoints;
      this->N = gridPoints.size(); // locations of particles in the domain
      // std::cout << "N: " << N << std::endl;
      // std::cout << "gridPoints.S: " << gridPoints.size() << std::endl;
      rhs = Vec(N);
      for (size_t i = 0; i < this->N; i++) {
        this->rhs(i) = F->RHSFunction(gridPoints[i]);
      }
    }
    #endif

      #ifdef noToLS
      int sqrtRootN;
      int L;
      void set_Uniform_Nodes() {
        std::vector<double> Nodes1D;
        for (int k=0; k<sqrtRootN; ++k) {
          Nodes1D.push_back(-L+2.0*L*(k+1.0)/(sqrtRootN+1.0));
        }
        pts2D temp1;
        for (int j=0; j<sqrtRootN; ++j) {
          for (int k=0; k<sqrtRootN; ++k) {
              temp1.x	=	Nodes1D[k];
              temp1.y	=	Nodes1D[j];
              this->gridPoints.push_back(temp1);
              particles_X.push_back(temp1);
              particles_Y.push_back(temp1);
          }
        }
      }

      void set_Standard_Cheb_Nodes() {
        std::vector<double> Nodes1D;
        for (int k=0; k<sqrtRootN; ++k) {
          Nodes1D.push_back(-cos((k+0.5)/sqrtRootN*PI));
        }
        pts2D temp1;
        for (int j=0; j<sqrtRootN; ++j) {
          for (int k=0; k<sqrtRootN; ++k) {
            temp1.x	=	Nodes1D[k];
            temp1.y	=	Nodes1D[j];
            gridPoints.push_back(temp1);
            particles_X.push_back(temp1);
            particles_Y.push_back(temp1);
          }
        }
      }

    userkernel(int sqrtRootN, double L, int nLevelsUniform) : kernel() {
      this->sqrtRootN = sqrtRootN;
      isTrans		=	true;
      isHomog		=	true;
      isLogHomog	=	false;
      alpha		=	-1.0;
      this->L = L;
      // this->set_Uniform_Nodes();
      this->set_Standard_Cheb_Nodes();
      this->N = gridPoints.size(); // locations of particles in the domain
      #ifdef USE_integralEqn
        double h = 1.0/sqrtRootN;
        // std::cout << "particles.size(): " << particles.size() << std::endl;
        double *a,*b;
        a = new double[2];
        b = new double[2];
        a[0] = 0;
        a[1] = 0;

        b[0] = h*0.5;
        b[1] = h*0.5;
        Kii = double_integral(a,b);
        Kii += 1; //For Second kind integral equation.
        h2 = 1.0/(sqrtRootN*sqrtRootN);
        #endif
  }
  #endif

  #ifdef USE_LS
  dtype getMatrixEntry(const unsigned i, const unsigned j) { //LS
    dtype output;
    output = F->getMatrixEntry(i, j);
    // return Aexplicit(i, j).real();
    return output;
    // return output.real()+I*output.imag();
  }
  #endif

  #ifdef USE_real
  dtype getMatrixEntry(const unsigned i, const unsigned j) {
    if (i==j)
      return pow(1000*this->N, 1.0/2.0);
      // return 1;
    else {
    	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
    	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
    	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
      double R = sqrt(R2);
      // return besselJ(0, R);
      // return 1.0/R;
      return 1.0/R;
    }
  }
  #endif

  #ifdef USE_integralEqn
  double Laplacian_2D(const unsigned i, const unsigned j) {
    pts2D r1 = particles_X[i];
    pts2D r2 = particles_X[j];
    double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
    double R = sqrt(R2);
    return log(R);
  }

  dtype getMatrixEntry(const unsigned i,const unsigned j){
  		// std::cout << "Kii: " << Kii << std::endl;
  		// std::cout << "h3: " << h3 << std::endl;
  		double res;
  		if(i==j)
  			res = Kii;//  triple_integral();
  		else{
  			res = h2 * Laplacian_2D(i,j);
  		}
  		return res;
  	}
  #endif

  #ifdef USE_Hankel
  dtype getMatrixEntry(const unsigned i, const unsigned j) {
    if (i==j)
      return pow(1000*this->N, 1.0/2.0);
      // return 1;
    else {
    	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
    	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
    	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
      double R = sqrt(R2);
      // return besselJ(0, R);
      // return 1.0/R;
      return I*(besselJ(0, R) + I*besselY(0, R))/4.0;
    }
  }
  #endif

  // dtype getMatrixEntry(const unsigned i, const unsigned j) {
  //   if (i==j)
  //     return pow(1000*this->N, 1.0/2.0);
  //     // return 1;
  //   else {
  //   	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
  //   	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
  //   	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
  //     double R = sqrt(R2);
  //     // return besselJ(0, R);
  //     // return 1.0/R;
  //     return I*(besselJ(0, R) + I*besselY(0, R))/4.0;
  //   }
  // }

  // dtype getMatrixEntry(const unsigned i, const unsigned j) {
  //   return(Aexplicit(i,j));
  //   // return(Aexplicit(i,j)+Aexplicit(j,i));
  // }


	// dtype getMatrixEntry(const unsigned i, const unsigned j) {
	// 	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
	// 	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
	// 	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
  //   double R = sqrt(R2);
  //   // return 1.0;
	// 	// if (R < 1e-10) {
	// 	// 	return 1.0;
	// 	// }
	// 	if (i==j) {
	// 		return 0.0;
	// 	}
	// 	else {
	// 		return 1.0/R;
	// 	}
	// }

  // dtype getMatrixEntry(const unsigned i, const unsigned j) {
	// 	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
	// 	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
	// 	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
	// 	if (R2 < 1e-10) {
	// 		return 1.0;
	// 	}
	// 	else if (R2 < a*a) {
	// 		return 0.5*R2*log(R2)/a/a;
	// 	}
	// 	else {
	// 		return 0.5*log(R2);
	// 	}
	// }

	// dtype getMatrixEntry(const unsigned i, const unsigned j) {
	// 	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
	// 	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
	// 	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
	// 	double R	=	sqrt(R2);
  //   // if (i==j) {
  //   //   return 1;//10.0*exp(I*1.0*R);
  //   // }
  //   	// if (R < a) {
  // 		// 	return R/a+1.0;
  // 		// }
  // 		// else {
  //     //   return exp(I*1.0*R)/R;
  // 		// }
  //     return exp(I*1.0*R);
	// }

  ~userkernel() {
    #ifdef yestoLS
      delete F;
    #endif
  };
};

// class userkernel: public kernel {
// private:
//   int nCones_LFR;//number of cones to be used in DAFMM
//   double TOL_POW;//currently not in use
//   int nLevelsUniform; //number of levels to be constructed in case of uniform tree
//   int yes2DFMM; //DFMM or FMM
//
// public:
//   // Mat Aexplicit;
//   LS2DTree* F;
//   std::vector<pts2D> particles_X;//dummy
//   std::vector<pts2D> particles_Y;//dummy
//   int N;
//   std::vector<pts2D> gridPoints; // location of particles in the domain
//   Vec rhs;
//   int nLevelsResolve;
//   userkernel(int nChebNodes, int treeAdaptivity, int degreeOfBases, double L, std::vector<pts2D>& particles_X, std::vector<pts2D>& particles_Y, int nLevelsUniform) : kernel(particles_X, particles_Y) {
//     isTrans		=	true;
//     isHomog		=	true;
//     isLogHomog	=	false;
//     alpha		=	-1.0;
//
//     this->nCones_LFR			=	16;//number of cones to be used in DAFMM
//     this->TOL_POW         = 9;//currently not in use
//     this->nLevelsUniform  = nLevelsUniform; //number of levels to be constructed in case of uniform tree
//     this->yes2DFMM				=	0; //DFMM or FMM
//     this->F	=	new LS2DTree(nCones_LFR, nChebNodes, L, yes2DFMM, TOL_POW, particles_X, particles_Y, kappa, degreeOfBases, treeAdaptivity, nLevelsUniform);
//     this->gridPoints = F->gridPoints;
//     this->nLevelsResolve = F->nLevels;
//     this->N = gridPoints.size(); // locations of particles in the domain
//     // std::cout << "N: " << N << std::endl;
//     // std::cout << "gridPoints.S: " << gridPoints.size() << std::endl;
//     rhs = Vec(N);
//     for (size_t i = 0; i < this->N; i++) {
//       this->rhs(i) = F->RHSFunction(gridPoints[i]);
//     }
//     // Aexplicit = Mat::Random(576,576);
//     // Mat B = Aexplicit.transpose();
//     // Aexplicit = Aexplicit+B;
//     // Aexplicit = Mat::Random(1008,1008);
//   }
//
//   // dtype getMatrixEntry(const unsigned i, const unsigned j) { //LS
//   //   dtype output;
//   //   output = F->getMatrixEntry(i, j);
//   //   // return Aexplicit(i, j).real();
//   //   return output;
//   //   // return output.real()+I*output.imag();
//   // }
//
//   dtype getMatrixEntry(const unsigned i, const unsigned j) {
//     if (i==j)
//       return pow(1000*this->N, 1.0/2.0);
//       // return 1;
//     else {
//     	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
//     	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
//     	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
//       double R = sqrt(R2);
//       // return besselJ(0, R);
//       // return 1.0/R;
//       return I*(besselJ(0, R) + I*besselY(0, R))/4.0;
//     }
//   }
//
//   // dtype getMatrixEntry(const unsigned i, const unsigned j) {
//   //   return(Aexplicit(i,j));
//   //   // return(Aexplicit(i,j)+Aexplicit(j,i));
//   // }
//
//
// 	// dtype getMatrixEntry(const unsigned i, const unsigned j) {
// 	// 	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
// 	// 	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
// 	// 	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
//   //   double R = sqrt(R2);
//   //   // return 1.0;
// 	// 	// if (R < 1e-10) {
// 	// 	// 	return 1.0;
// 	// 	// }
// 	// 	if (i==j) {
// 	// 		return 0.0;
// 	// 	}
// 	// 	else {
// 	// 		return 1.0/R;
// 	// 	}
// 	// }
//
//   // dtype getMatrixEntry(const unsigned i, const unsigned j) {
// 	// 	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
// 	// 	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
// 	// 	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
// 	// 	if (R2 < 1e-10) {
// 	// 		return 1.0;
// 	// 	}
// 	// 	else if (R2 < a*a) {
// 	// 		return 0.5*R2*log(R2)/a/a;
// 	// 	}
// 	// 	else {
// 	// 		return 0.5*log(R2);
// 	// 	}
// 	// }
//
// 	// dtype getMatrixEntry(const unsigned i, const unsigned j) {
// 	// 	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
// 	// 	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
// 	// 	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
// 	// 	double R	=	sqrt(R2);
//   //   // if (i==j) {
//   //   //   return 1;//10.0*exp(I*1.0*R);
//   //   // }
//   //   	// if (R < a) {
//   // 		// 	return R/a+1.0;
//   // 		// }
//   // 		// else {
//   //     //   return exp(I*1.0*R)/R;
//   // 		// }
//   //     return exp(I*1.0*R);
// 	// }
//
//   ~userkernel() {
//     delete F;
//   };
// };

// class userkernel: public kernel {
// public:
// 	dtype RHSFunction(const pts2D r) {
// 		dtype q = -1.0+I;
// 		return q;
// 	};
// 	userkernel(std::vector<pts2D>& particles_X, std::vector<pts2D>& particles_Y): kernel(particles_X, particles_Y) {
// 		isTrans		=	true;
// 		isHomog		=	true;
// 		isLogHomog	=	false;
// 		alpha		=	-1.0;
// 	};
// 	// dtype getMatrixEntry(const unsigned i, const unsigned j) {
// 	// 	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
// 	// 	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
// 	// 	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
// 	// 	double R	=	sqrt(R2);
// 	// 	if (R < a) {
// 	// 		return R/a+1.0;
// 	// 	}
// 	// 	else {
// 	// 		return (1.0+I)*a/R;
// 	// 	}
// 	// }
//
// 	dtype getMatrixEntry(const unsigned i, const unsigned j) {
// 		pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
// 		pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
// 		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
// 		double R	=	sqrt(R2);
// 		return exp(I*1.0*R);
// 	}
//
// 	// #elif LOGR
// 	// userkernel(std::vector<pts2D> particles_X, std::vector<pts2D> particles_Y): kernel(particles_X, particles_Y) {
// 	// 	isTrans		=	true;
// 	// 	isHomog		=	false;
// 	// 	isLogHomog	=	true;
// 	// 	alpha		=	1.0;
// 	// };
// 	// double getMatrixEntry(const unsigned i, const unsigned j) {
// 	// 	pts2D r1 = particles_X[i];//particles_X is a member of base class FMM_Matrix
// 	// 	pts2D r2 = particles_X[j];//particles_X is a member of base class FMM_Matrix
// 	// 	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
// 	// 	if (R2 < 1e-10) {
// 	// 		return 1.0;
// 	// 	}
// 	// 	else if (R2 < a*a) {
// 	// 		return 0.5*R2*log(R2)/a/a;
// 	// 	}
// 	// 	else {
// 	// 		return 0.5*log(R2);
// 	// 	}
// 	// }
// 	// #endif
// 	~userkernel() {};
// };
#endif
