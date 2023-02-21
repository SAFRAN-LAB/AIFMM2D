// typedef std::complex<double> dtype;
class classGMRES {
public:
  classGMRES(){};
/// Calculate the Given rotation matrix----%%
 void givensRotation(dtype v1, dtype v2, dtype& cs, dtype& sn) {
   Vec temp(2);
   temp(0) = v1;
   temp(1) = v2;
   dtype t = temp.norm();
   cs = v1 / t;
   sn = v2 / t;
 }

void applyGivensRotation(Vec& h, Vec cs, Vec sn, int k, dtype& cs_k, dtype& sn_k) {
  // apply for ith column
  for (int i = 0; i < k; i++) {
    dtype temp = cs(i) * h(i) + sn(i) * h(i + 1);
    h(i+1) = -sn(i) * h(i) + cs(i) * h(i + 1);
    h(i)   = temp;
  }
  // update the next sin cos values for rotation
  givensRotation(h(k), h(k+1), cs_k, sn_k);
  // eliminate H(i + 1, i)
  h(k) = cs_k * h(k) + sn_k * h(k+1);
  h(k+1) = 0.0;
}

template <typename MatVecAccelerator>
void arnoldi(MatVecAccelerator* M, Mat Q, int k, Vec &h, Vec &q) {
  Vec temp = Q.col(k);
  M->MatVecProduct(temp, q); // Krylov Vector
  for (size_t i = 0; i < k+1; i++) { // Modified Gram-Schmidt, keeping the Hessenberg matrix
    h(i) = q.dot(Q.col(i));
    q = q - h(i) * Q.col(i);
  }
  h(k+1) = q.norm();
  q = q / h(k+1);
}

template <typename MatVecAccelerator, typename Preconditioner>
void arnoldi(MatVecAccelerator* M, Preconditioner* P , Mat Q, int k, Vec &h, Vec &q) {
  Vec temp;
  Vec temp1 = Q.col(k);
  M->MatVecProduct(temp1, temp);
  q = P->solve(temp); // Krylov Vector
  for (size_t i = 0; i < k+1; i++) { // Modified Gram-Schmidt, keeping the Hessenberg matrix
    h(i) = q.dot(Q.col(i));
    q = q - h(i) * Q.col(i);
  }
  h(k+1) = q.norm();
  q = q / h(k+1);
}

template <typename MatVecAccelerator>
void gmres(MatVecAccelerator* M, Vec b, int maxIterations, double threshold, Vec& x, double& error, int& noOfIterations, std::vector<double>& e) {
  int N = b.size(); //assuming A is square matrix
  int m = maxIterations;

  // use x as the initial vector
  x = Vec::Zero(N);
  Vec temp;
  M->MatVecProduct(x, temp);
  Vec r = b - temp;
  double b_norm = b.norm();
  error = r.norm()/b.norm();

  // initialize the 1D vectors
  Vec sn = Vec::Zero(m);
  Vec cs = Vec::Zero(m);
  Vec e1 = Vec::Zero(m+1);
  e1(0) = 1;
  // std::vector<double> e;
  e.push_back(error);
  double r_norm = r.norm();
  Mat Q = Mat::Zero(N,m+1);
  Mat H = Mat::Zero(m+1,m);
  Q.col(0) = r/r_norm;
  Vec beta = r_norm*e1;
  int K=m-1;

  for (size_t k = 0; k < m; k++) {
    // std::cout << "k: " << k << std::endl;
    // run arnoldi
    Vec h(k+2);
    Vec q;
    arnoldi(M, Q.block(0,0,N,k+1), k, h, q);
    H.block(0,k,k+2,1) = h;
    Q.block(0,k+1,N,1) = q;
    // eliminate the last element in H ith row and update the rotation matrix
    dtype cs_k, sn_k;
    applyGivensRotation(h, cs, sn, k, cs_k, sn_k);
    H.block(0,k,k+2,1) = h;
    cs(k) = cs_k;
    sn(k) = sn_k;
    // update the residual vector
    beta(k+1) = -sn(k) * beta(k);
    beta(k)     = cs(k) * beta(k);
    error       = abs(beta(k+1)) / b_norm;
    // std::cout << "k: " << k <<  " error: " << error << std::endl;
    // save the error
    e.push_back(error);
    if (error <= threshold) {
      // std::cout << "threshold: " << threshold << std::endl;
      // std::cout << "error: " << error << std::endl;
      K = k;
      break;
    }
  }
  noOfIterations = K+1;
  // if threshold is not reached, k = m at this point (and not m+1)
  // std::cout << "No. of iterations: " << K+1 << std::endl;
  // calculate the result
  Mat H_final = H.block(0,0,K+1,K+1);
  Vec y = H_final.colPivHouseholderQr().solve(beta.segment(0,K+1));
  x = x + Q.block(0,0,N,K+1) * y;
  // std::cout << "No. of iterations: " << K+1 << std::endl;
  // Vec errVec = A*x-b;
  // double errGMRES = errVec.norm()/b.norm();
  // std::cout << "Err: " << errGMRES << std::endl;
}

template <typename MatVecAccelerator, typename Preconditioner>
void gmres(MatVecAccelerator* M, Preconditioner* P, Vec b_old, int maxIterations, double threshold, Vec& x, double& error, int& noOfIterations, std::vector<double>& e) {
  int N = b_old.size(); //assuming A is square matrix
  int m = maxIterations;
  Vec b = P->solve(b_old);
  // use x as the initial vector
  x = Vec::Zero(N);
  Vec temp;
  M->MatVecProduct(x, temp);
  Vec r = b - P->solve(temp);
  double b_norm = b.norm();
  error = r.norm()/b.norm();

  // initialize the 1D vectors
  Vec sn = Vec::Zero(m);
  Vec cs = Vec::Zero(m);
  Vec e1 = Vec::Zero(m+1);
  e1(0) = 1;
  // std::vector<double> e;
  e.push_back(error);
  double r_norm = r.norm();
  Mat Q = Mat::Zero(N,m+1);
  Mat H = Mat::Zero(m+1,m);
  Q.col(0) = r/r_norm;
  Vec beta = r_norm*e1;
  int K=m-1;

  for (size_t k = 0; k < m; k++) {
    // std::cout << "k: " << k << std::endl;
    // run arnoldi
    Vec h(k+2);
    Vec q;
    arnoldi(M, P, Q.block(0,0,N,k+1), k, h, q);
    H.block(0,k,k+2,1) = h;
    Q.block(0,k+1,N,1) = q;
    // eliminate the last element in H ith row and update the rotation matrix
    dtype cs_k, sn_k;
    applyGivensRotation(h, cs, sn, k, cs_k, sn_k);
    H.block(0,k,k+2,1) = h;
    cs(k) = cs_k;
    sn(k) = sn_k;
    // update the residual vector
    beta(k+1) = -sn(k) * beta(k);
    beta(k)     = cs(k) * beta(k);
    error       = abs(beta(k+1)) / b_norm;
    // std::cout << "k: " << k <<  " error: " << error << std::endl;
    // save the error
    e.push_back(error);
    if (error <= threshold) {
      // std::cout << "threshold: " << threshold << std::endl;
      // std::cout << "error: " << error << std::endl;
      K = k;
      break;
    }
  }
  noOfIterations = K+1;
  // if threshold is not reached, k = m at this point (and not m+1)
  // std::cout << "No. of iterations: " << K+1 << std::endl;
  // calculate the result
  Mat H_final = H.block(0,0,K+1,K+1);
  Vec y = H_final.colPivHouseholderQr().solve(beta.segment(0,K+1));
  x = x + Q.block(0,0,N,K+1) * y;
  // std::cout << "No. of iterations: " << K+1 << std::endl;
  // Vec errVec = A*x-b;
  // double errGMRES = errVec.norm()/b.norm();
  // std::cout << "Err: " << errGMRES << std::endl;
}

~classGMRES(){};

};
