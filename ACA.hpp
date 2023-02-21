//
//  ACA.hpp
//
//
//  Created by Vaishnavi Gujjula on 12/4/20.
//
//
class LowRank {
public:
	userkernel* K;
	double tol_ACA;
	double tol_robust;
	double tolerance_or_rank;
	std::vector<int> row_indices;
	std::vector<int> col_indices;

	LowRank(userkernel* K, double tolerance_or_rank, std::vector<int>& row_indices, std::vector<int>& col_indices) {
		this->K					=	K;
		this->tolerance_or_rank	=	tolerance_or_rank;
		this->tol_ACA	=	tolerance_or_rank; //assuming it is tolerance
		// this->tol_ACA = pow(10,-1.0*tolerance_or_rank);
		this->row_indices = row_indices;
		this->col_indices = col_indices;
		this->tol_robust = pow(10.0,-10.0); //best is -10 for 16384 for gaussian
		//best is -7 for 36864 for Lens
	}

	void maxAbsVector(const Vec& v, const std::set<int>& allowed_indices,
																dtype& max, int& index
															 ) {
		std::set<int>::iterator it;
		index = *allowed_indices.begin();
		max   = v(index);

		for(it = allowed_indices.begin(); it != allowed_indices.end(); it++) {
				if(abs(v(*it)) > abs(max)) {
						index   =   *it;
						max     =   v(index);
				}
		}
	}

	void rookPiv(std::vector<int>& row_ind, std::vector<int>& col_ind, int& computed_rank, Mat& L, Mat& R)
	{
	    // Indices which have been used:
	    int n_rows = row_indices.size();
			int n_cols = col_indices.size();
			// std::cout << "n_rows: " << n_rows << "	n_cols: " << n_cols << std::endl;
			if (n_rows == 0 || n_cols == 0) {
				computed_rank = 0;
				return;
			}
	    // Indices that are remaining:
	    std::set<int> remaining_row_ind;
	    std::set<int> remaining_col_ind;

	    // Bases:
	    std::vector<Vec> u;
	    std::vector<Vec> v;
			std::vector<Vec> u_original;
	    std::vector<Vec> v_original;

	    for(int k = 0; k < n_rows; k++)
	    {
	        remaining_row_ind.insert(k);
	    }

	    for(int k = 0; k < n_cols; k++)
	    {
	        remaining_col_ind.insert(k);
	    }

	    dtype max, gamma;

	    // Initialize the matrix norm and the the first row index
	    double matrix_norm = 0;
	    row_ind.push_back(0);
	    remaining_row_ind.erase(0);

	    // Stores the pivot entry of the considered row / col:
	    int pivot;

	    int target_rank = 0;
	    // This would get updated:
	    computed_rank = 0;
	    Vec row, col;
			Vec row_temp, col_temp;

	    double tolerance = 0;
	    // These quantities in finding the stopping criteria:
	    double row_squared_norm, row_norm, col_squared_norm, col_norm;

	    if(tolerance_or_rank < 1)
	        tolerance = tolerance_or_rank;
	    else
	        target_rank = tolerance_or_rank;

	    // So these would be particularly useful for poorly conditioned matrices:
	    int max_tries = 10;
	    int count;

	    // Repeat till the desired tolerance / rank is obtained
	    do
	    {
					// std::cout << "in while" << std::endl;
	        // Generation of the row
	        // Row of the residuum and the pivot column
	        // By calling row_ind.back(), we are getting the last pushed number
	        row = K->getRow(row_indices[row_ind.back()], col_indices);
					row_temp = row;
	        for(int i = 0; i < computed_rank; i++)
	        {
	            row = row - u[i](row_ind.back()) * v[i];
	        }

	        this->maxAbsVector(row, remaining_col_ind, max, pivot);
	        count = 0;

	        // Alternating upon each call:
	        bool eval_at_end = false;
	        // Toggling randomness
	        bool use_randomization = true;

	        // This while loop is needed if in the middle of the algorithm the
	        // row happens to be exactly the linear combination of the previous rows
	        // upto some tolerance. i.e. prevents from ACA throwing false positives
	        while (fabs(max) < tolerance &&
	               count < max_tries   &&
	               remaining_col_ind.size() > 0 &&
	               remaining_row_ind.size() > 0
	              )
	        {
	            row_ind.pop_back();
	            int new_row_ind;

	            // When rank < 3, we will just choose entries from the ends of the matrix:
	            if(computed_rank < 3)
	            {
	                if(eval_at_end == true)
	                {
	                    new_row_ind = *remaining_row_ind.end();
	                }

	                else
	                {
	                    new_row_ind = *remaining_row_ind.begin();
	                }

	                eval_at_end = !(eval_at_end);
	            }

	            // However, when we have rank >=3, we will choose the entries such that
	            // the newly picked entry is at the mid-point of the already chosen ones:
	            else
	            {
	                if(use_randomization == true)
	                {
	                    std::set<int>::const_iterator it(remaining_row_ind.begin());
	                    std::advance(it, rand() % remaining_row_ind.size());
	                    new_row_ind = *it;
	                }

	                else
	                {
	                    std::vector<int> row_ind_sort(row_ind);
	                    std::sort(row_ind_sort.begin(), row_ind_sort.end());
	                    std::vector<int> row_ind_diff(row_ind_sort.size() - 1);

	                    int max = 0;
	                    int idx = 0;

	                    for(int i = 0; i < row_ind_sort.size() - 1; i++)
	                    {
	                        row_ind_diff[i] = row_ind_sort[i+1] - row_ind_sort[i];
	                        if(row_ind_diff[i] > max)
	                        {
	                            idx = i;
	                            max = row_ind_diff[i];
	                        }
	                    }

	                    new_row_ind = row_ind_sort[idx] + max / 2;
	                }

	                use_randomization = !(use_randomization);
	            }

	            row_ind.push_back(new_row_ind);
	            remaining_row_ind.erase(new_row_ind);
	            // Generation of the row
	            // Row of the residuum and the pivot column
							row = K->getRow(row_indices[new_row_ind], col_indices);
							row_temp = row;
							for(int i = 0; i < computed_rank; i++)
	            {
	                row = row - u[i](row_ind.back()) * v[i];
	            }

	            this->maxAbsVector(row, remaining_col_ind, max, pivot);
	            count++;
	        }

	        // In case it failed to resolve in the previous step,
	        // we break out of the dowhile loop:
	        if (count == max_tries ||
	            remaining_col_ind.size() == 0 ||
	            remaining_row_ind.size() == 0
	           )
	        {
	            break;
	        }

	        // Resetting count back to zero for columns:
	        count = 0;

	        col_ind.push_back(pivot);
	        remaining_col_ind.erase(pivot);
	        // Normalizing constant
	        gamma = double(1.0) / max;

	        // Generation of the column
	        // Column of the residuum and the pivot row
					col = K->getCol(col_indices[col_ind.back()], row_indices);
					col_temp = col;
					for(int i = 0; i < computed_rank; i++)
	        {
	            col = col - v[i](col_ind.back()) * u[i];
	        }

	        this->maxAbsVector(col, remaining_row_ind, max, pivot);
	        // Repeating the same randomization we carried out for the rows, now for the columns:
	        while (fabs(max)<tolerance &&
	               count < max_tries &&
	               remaining_col_ind.size() >0 &&
	               remaining_row_ind.size() >0
	              )
	        {
	            col_ind.pop_back();

	            int new_col_ind;

	            if(col_ind.size() < 3)
	            {
	                if(eval_at_end)
	                {
	                    new_col_ind = *remaining_col_ind.end();
	                }

	                else
	                {
	                    new_col_ind = *remaining_col_ind.begin();
	                }

	                eval_at_end = !eval_at_end;
	            }

	            else
	            {
	                if(use_randomization == true)
	                {
	                    std::set<int>::const_iterator it(remaining_col_ind.begin());
	                    std::advance(it, rand() % remaining_col_ind.size());
	                    new_col_ind = *it;
	                }

	                else
	                {
	                    std::vector<int> col_ind_sort(col_ind);
	                    std::sort(col_ind_sort.begin(), col_ind_sort.end());
	                    std::vector<int> col_ind_diff(col_ind_sort.size() - 1);

	                    int max = 0;
	                    int idx = 0;

	                    for(int i = 0; i < col_ind_sort.size() - 1; i++)
	                    {
	                        col_ind_diff[i] = col_ind_sort[i+1] - col_ind_sort[i];
	                        if(col_ind_diff[i] > max)
	                        {
	                            idx = i;
	                            max = col_ind_diff[i];
	                        }
	                    }

	                    new_col_ind = col_ind_sort[idx] + max / 2;
	                }

	                use_randomization = !(use_randomization);
	            }

	            col_ind.push_back(new_col_ind);
	            remaining_col_ind.erase(new_col_ind);

	            // Generation of the column
	            // Column of the residuum and the pivot row:
							col = K->getCol(col_indices[new_col_ind], row_indices);
							col_temp = col;
							for(int i = 0; i < computed_rank; i++)
	            {
	                col = col - v[i](col_ind.back()) * u[i];
	            }

	            this->maxAbsVector(col, remaining_row_ind, max, pivot);
	            count++;
	        }

	        row_ind.push_back(pivot);
	        remaining_row_ind.erase(pivot);

	        // New vectors
	        u.push_back(gamma * col);
	        v.push_back(row);
					u_original.push_back(col_temp);
					v_original.push_back(row_temp);

	        // New approximation of matrix norm
	        row_squared_norm = row.squaredNorm();
	        row_norm         = sqrt(row_squared_norm);

	        col_squared_norm = col.squaredNorm();
	        col_norm         = sqrt(col_squared_norm);

	        // Updating the matrix norm:
	        matrix_norm += std::abs(gamma * gamma * row_squared_norm * col_squared_norm);

	        for(int j = 0; j < computed_rank; j++)
	        {
	            matrix_norm += 2.0 * std::abs(u[j].dot(u.back()))
	                               * std::abs(v[j].dot(v.back()));
	        }

	        computed_rank++;
	    }
	    while(((tolerance_or_rank < 1) ?
	          computed_rank * (n_rows + n_cols) * row_norm * col_norm >
	          fabs(max) * tolerance * matrix_norm
	          : computed_rank < target_rank)
	          &&
	          computed_rank < fmin(n_rows, n_cols)+1
	         );
		  // std::cout << "out while" << std::endl;
		 	// If the computed_rank is >= to full-rank
	    // then return the trivial full-rank decomposition
			row_ind.pop_back(); // removing the last pivot when while has ended
			if (computed_rank >= fmin(n_rows, n_cols))
	    {
					// std::cout << "in full rank; " << u_original.size() << ", " << v_original.size() << std::endl;
					if (n_rows < n_cols)
	        {
						computed_rank = n_rows;
	        }
	        else
	        {
	            computed_rank = n_cols;
	        }
	    }
			// std::cout << "computed_rank: " << computed_rank << "	row_ind: " << row_ind.size() << "	col_ind: " << col_ind.size() << std::endl;
	    // This is when ACA has succeeded:
      L = Mat(n_rows, computed_rank);
      R = Mat(computed_rank, n_cols);

      for (int j = 0; j < computed_rank; j++)
      {
          L.col(j) = u_original[j];
          R.row(j) = v_original[j];
      }
	}

	void ACA_only_nodes(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, Mat &Ac, Mat &Ar, Mat& L, Mat& R) {
		int col_index;
		int row_index;
		int N1 = row_indices.size();
		int N2 = col_indices.size();
		Vec row(N2), col(N1), v(N2), u(N1);
		Vec row_temp, col_temp;
		std::vector<Vec> Uvec;
		std::vector<Vec> Vvec;
		std::vector<Vec> AcVec;
		std::vector<Vec> ArVec;
		computed_rank = 0;
		dtype max;
		int min = N1;
		if (N1 > N2) {
			min = N2;
		}
		// Mat Ac = Mat(N1, min);
		// Mat Ar = Mat(min, N2);

		std::set<int> remaining_row_ind;
		std::set<int> remaining_col_ind;

		for(int l = 0; l < N1; l++) {
				remaining_row_ind.insert(l);
		}
		for(int l= 0; l < N2; l++) {
				remaining_col_ind.insert(l);
		}
		if (N1 < N2) {
			int l_local;
			for(l_local = 0; l_local < N1; l_local++) {
				row_index=l_local;
				row = K->getRow(row_indices[row_index], col_indices);
				this->maxAbsVector(row, remaining_col_ind, max, col_index);
				if(abs(row(int(col_index))) > pow(10,-16)) {
					break;
				}
			}
			if (l_local == N1) {
				Ac = Mat(N1,computed_rank);
				Ar = Mat(computed_rank,N2);
				// Ac_modified = Ac.block(0,0,N1,computed_rank);
				// Ar_modified = Ar.block(0,0,computed_rank,N2);
				return;
			}
			v=row;
			row_bases.push_back(row_index);
			col_bases.push_back(int(col_index));
			col = K->getCol(col_indices[col_index], row_indices);
			u	=	col/row(int(col_index));
			Uvec.push_back(u);
			Vvec.push_back(v);
			// Ac.col(computed_rank) = col;
			// Ar.row(computed_rank) = row;
			AcVec.push_back(col);
			ArVec.push_back(row);
			remaining_col_ind.erase(col_index);
			remaining_row_ind.erase(row_index);
			computed_rank = 1;

			double normS	=	0.0;
			double prev_normS = DBL_MAX;
			this->maxAbsVector(col, remaining_row_ind, max, row_index);

			while (abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
			// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && (abs(prev_normS-normS) >= tol_ACA*prev_normS || abs(prev_normS-normS) >= tol_ACA) && computed_rank < min) {
			// while(u.norm()*v.norm() > tol_ACA*normS &&  computed_rank < min) {
			// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
				row_bases.push_back(int(row_index));
				row = K->getRow(row_indices[row_index], col_indices);
				row_temp = row;
				for (int l=0; l<=computed_rank-1; ++l){
					for (size_t s = 0; s < Vvec[l].size(); s++) {
						row(s)	-=	Uvec[l](row_index)*Vvec[l](s);
					}
				}
				this->maxAbsVector(row, remaining_col_ind, max, col_index);
				col_bases.push_back(int(col_index));
				v = row;
				col = K->getCol(col_indices[col_index], row_indices);
				col_temp = col;
				for (int l=0; l<=computed_rank-1; ++l){
					for (size_t s = 0; s < Uvec[l].size(); s++) {
						col(s)	-=	Vvec[l](col_index)*Uvec[l](s);
					}
				}
				u	=	col/row(int(col_index));
				// if(u.norm()< tol_robust || v.norm()< tol_robust) {
				// 	row_bases.pop_back();
				// 	col_bases.pop_back();
				// 	// std::cout << "break" << std::endl;
				// 	break;
				// }
				Uvec.push_back(u);
				Vvec.push_back(v);
				// Ac.col(computed_rank) = col_temp;
				// Ar.row(computed_rank) = row_temp;
				AcVec.push_back(col_temp);
				ArVec.push_back(row_temp);

				++computed_rank;
				remaining_col_ind.erase(col_index);
				remaining_row_ind.erase(row_index);
				if (computed_rank != 2)
					prev_normS = normS;
				normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
				for (int l=0; l<=computed_rank-1; ++l){
					normS	+=	2*abs(Uvec[l].dot(u) * Vvec[l].dot(v));
				}
				normS	=	sqrt(normS);
				this->maxAbsVector(col, remaining_row_ind, max, row_index);
				// std::cout << "check: " << (u.norm()*v.norm() > tol_ACA*normS) << std::endl;
				// std::cout << "uv: " << u.norm()*v.norm() << std::endl;
				// std::cout << "ts: " << tol_ACA*normS << std::endl;
				// std::cout << "tol_ACA: " << tol_ACA << std::endl;
			}
			//////////
			// for (size_t i = 0; i < Uvec.size(); i++) {
			// 	std::cout << "i: " << i << std::endl << Uvec[i] << std::endl;
			// }
			// for (size_t i = 0; i < Vvec.size(); i++) {
			// 	std::cout << "i: " << i << std::endl << Vvec[i] << std::endl;
			// }
			L = Mat::Zero(computed_rank, computed_rank);
			R = Mat::Zero(computed_rank, computed_rank);
			if(computed_rank > 0) {
				for (size_t i = 0; i < computed_rank; i++) {
					L(i,i) = 1.0;
					if (i >= 1) {
						for (size_t j = 0; j <= i-1; j++) {
							L(i,j) = Uvec[j](row_bases[i]);//Vvec[j](col_bases[i]);
						}
					}
				}
				for (size_t i = 0; i < computed_rank; i++) {
					R(i,i) = Vvec[i](col_bases[i]);//Uvec[i](row_bases[i]);
					if (i >= 1) {
						for (size_t j = 0; j <= i-1; j++) {
							R(j,i) = Vvec[j](col_bases[i]);//Uvec[i](row_bases[i]);
						}
					}
				}
			 }

		}
		else {
			int l_local;
			for(l_local = 0; l_local < N2; l_local++) {
				col_index=l_local;
				col = K->getCol(col_indices[col_index], row_indices);
				this->maxAbsVector(col, remaining_row_ind, max, row_index);
				if(abs(col(int(row_index))) > pow(10,-16)) {
					break;
				}
			}
			if (l_local == N2) {
				// Ac_modified = Ac.block(0,0,N1,computed_rank);
				// Ar_modified = Ar.block(0,0,computed_rank,N2);
				Ac = Mat(N1,computed_rank);
				Ar = Mat(computed_rank,N2);
				return;
			}

			col_bases.push_back(col_index);
			v=col;
			row_bases.push_back(int(row_index));
			row = K->getRow(row_indices[row_index], col_indices);
			u	=	row/col(int(row_index));
			Uvec.push_back(u);
			Vvec.push_back(v);
			// Ac.col(computed_rank) = col;
			// Ar.row(computed_rank) = row;
			AcVec.push_back(col);
			ArVec.push_back(row);
			computed_rank = 1;

			remaining_row_ind.erase(row_index);
			remaining_col_ind.erase(col_index);

			double normS	=	0.0;
			double prev_normS = DBL_MAX;
			this->maxAbsVector(row, remaining_col_ind, max, col_index);

			// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && (abs(prev_normS-normS) >= tol_ACA*prev_normS || abs(prev_normS-normS) >= tol_ACA) && computed_rank < min) {
			// while(u.norm()*v.norm() > tol_ACA*normS &&  computed_rank < min) {
			// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
			while (abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
				col_bases.push_back(int(col_index));
			  col = K->getCol(col_indices[col_index], row_indices);
				col_temp = col;
			  for (int l=0; l<=computed_rank-1; ++l){
			    for (size_t s = 0; s < Vvec[l].size(); s++) {
			      col(s)	-=	Uvec[l](col_index)*Vvec[l](s);
			    }
			  }
			  this->maxAbsVector(col, remaining_row_ind, max, row_index);
			  row_bases.push_back(int(row_index));
			  v=col;
			  row = K->getRow(row_indices[row_index], col_indices);
				row_temp = row;
			  for (int l=0; l<=computed_rank-1; ++l){
			    for (size_t s = 0; s < Uvec[l].size(); s++) {
			      row(s)	-=	Vvec[l](row_index)*Uvec[l](s);
			    }
			  }
			  u	=	row/col(int(row_index));
			  // if(u.norm()< tol_robust || v.norm()< tol_robust) {
			  //   col_bases.pop_back();
			  //   row_bases.pop_back();
			  //   break;
			  // }
			  Uvec.push_back(u);
			  Vvec.push_back(v);
				// Ac.col(computed_rank) = col_temp;
				// Ar.row(computed_rank) = row_temp;
				AcVec.push_back(col_temp);
				ArVec.push_back(row_temp);

			  ++computed_rank;
				remaining_row_ind.erase(row_index);
			  remaining_col_ind.erase(col_index);
			  if (computed_rank != 2)
			    prev_normS = normS;
			  normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
			  for (int l=0; l<=computed_rank-1; ++l){
			    normS	+=	2*abs(Uvec[l].dot(u) * Vvec[l].dot(v));
			  }
			  normS	=	sqrt(normS);
			  this->maxAbsVector(row, remaining_col_ind, max, col_index);
			}//while
			//////////
			// for (size_t i = 0; i < Uvec.size(); i++) {
			// 	std::cout << "i: " << i << std::endl << Uvec[i] << std::endl;
			// }
			// for (size_t i = 0; i < Vvec.size(); i++) {
			// 	std::cout << "i: " << i << std::endl << Vvec[i] << std::endl;
			// }
			L = Mat::Zero(computed_rank, computed_rank);
			R = Mat::Zero(computed_rank, computed_rank);
			if(computed_rank > 0) {
				for (size_t i = 0; i < computed_rank; i++) {
					L(i,i) = 1.0;
					if (i >= 1) {
						for (size_t j = 0; j <= i-1; j++) {
							L(j,i) = Uvec[j](col_bases[i]);//Vvec[j](row_bases[i]);
						}
					}
				}
				for (size_t i = 0; i < computed_rank; i++) {
					R(i,i) = Vvec[i](row_bases[i]);//Uvec[i](col_bases[i]);
					if (i >= 1) {
						for (size_t j = 0; j <= i-1; j++) {
							R(i,j) = Vvec[j](row_bases[i]);//Uvec[j](col_bases[i]);
						}
					}
				}
			 }
		}//else

		// Ac_modified = Ac.block(0,0,N1,computed_rank);
		// Ar_modified = Ar.block(0,0,computed_rank,N2);
		Ac = Mat(N1,computed_rank);
		Ar = Mat(computed_rank,N2);
		for (size_t i = 0; i < computed_rank; i++) {
			Ac.col(i) = AcVec[i];
			Ar.row(i) = ArVec[i];
		}
	}

	void ACA_only_nodes(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, Mat &Ac, Mat &Ar) {
		int col_index;
		int row_index;
		int N1 = row_indices.size();
		int N2 = col_indices.size();
		Vec row(N2), col(N1), v(N2), u(N1);
		Vec row_temp, col_temp;
		std::vector<Vec> Uvec;
		std::vector<Vec> Vvec;
		std::vector<Vec> AcVec;
		std::vector<Vec> ArVec;
		computed_rank = 0;
		dtype max;
		int min = N1;
		if (N1 > N2) {
			min = N2;
		}
		// Mat Ac = Mat(N1, min);
		// Mat Ar = Mat(min, N2);

		std::set<int> remaining_row_ind;
		std::set<int> remaining_col_ind;

		for(int l = 0; l < N1; l++) {
				remaining_row_ind.insert(l);
		}
		for(int l= 0; l < N2; l++) {
				remaining_col_ind.insert(l);
		}
		if (N1 < N2) {
			int l_local;
			for(l_local = 0; l_local < N1; l_local++) {
				row_index=l_local;
				row = K->getRow(row_indices[row_index], col_indices);
				this->maxAbsVector(row, remaining_col_ind, max, col_index);
				if(abs(row(int(col_index))) > pow(10, -16)) {
					break;
				}
			}
			if (l_local == N1) {
				Ac = Mat(N1,computed_rank);
				Ar = Mat(computed_rank,N2);
				// Ac_modified = Ac.block(0,0,N1,computed_rank);
				// Ar_modified = Ar.block(0,0,computed_rank,N2);
				return;
			}
			v=row/row(int(col_index));
			row_bases.push_back(row_index);
			col_bases.push_back(int(col_index));
			col = K->getCol(col_indices[col_index], row_indices);
			u	=	col;
			Uvec.push_back(u);
			Vvec.push_back(v);
			// Ac.col(computed_rank) = col;
			// Ar.row(computed_rank) = row;
			AcVec.push_back(col);
			ArVec.push_back(row);
			remaining_col_ind.erase(col_index);
			remaining_row_ind.erase(row_index);
			computed_rank = 1;

			double normS	=	0.0;
			double prev_normS = DBL_MAX;
			this->maxAbsVector(col, remaining_row_ind, max, row_index);

			while (abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
			// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && (abs(prev_normS-normS) >= tol_ACA*prev_normS || abs(prev_normS-normS) >= tol_ACA) && computed_rank < min) {
			// while(u.norm()*v.norm() > tol_ACA*normS &&  computed_rank < min) {
			// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
				row_bases.push_back(int(row_index));
				row = K->getRow(row_indices[row_index], col_indices);
				row_temp = row;
				for (int l=0; l<=computed_rank-1; ++l){
					for (size_t s = 0; s < Vvec[l].size(); s++) {
						row(s)	-=	Uvec[l](row_index)*Vvec[l](s);
					}
				}
				this->maxAbsVector(row, remaining_col_ind, max, col_index);
				col_bases.push_back(int(col_index));
				for (size_t l = 0; l < row.size(); l++) {
					v(l) = row(l)/row(int(col_index));
				}
				col = K->getCol(col_indices[col_index], row_indices);
				col_temp = col;
				for (int l=0; l<=computed_rank-1; ++l){
					for (size_t s = 0; s < Uvec[l].size(); s++) {
						col(s)	-=	Vvec[l](col_index)*Uvec[l](s);
					}
				}
				u	=	col;
				// if(u.norm()< tol_robust || v.norm()< tol_robust) {
				// 	row_bases.pop_back();
				// 	col_bases.pop_back();
				// 	break;
				// }
				Uvec.push_back(u);
				Vvec.push_back(v);
				// Ac.col(computed_rank) = col_temp;
				// Ar.row(computed_rank) = row_temp;
				AcVec.push_back(col_temp);
				ArVec.push_back(row_temp);

				++computed_rank;
				remaining_col_ind.erase(col_index);
				remaining_row_ind.erase(row_index);
				if (computed_rank != 2)
					prev_normS = normS;
				normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
				for (int l=0; l<=computed_rank-1; ++l){
					normS	+=	2*abs(Uvec[l].dot(u) * Vvec[l].dot(v));
				}
				normS	=	sqrt(normS);
				this->maxAbsVector(col, remaining_row_ind, max, row_index);
			}
		}
		else {
			int l_local;
			for(l_local = 0; l_local < N2; l_local++) {
				col_index=l_local;
				col = K->getCol(col_indices[col_index], row_indices);
				this->maxAbsVector(col, remaining_row_ind, max, row_index);
				if(abs(col(int(row_index))) > pow(10,-16)) {
					break;
				}
			}
			if (l_local == N2) {
				// Ac_modified = Ac.block(0,0,N1,computed_rank);
				// Ar_modified = Ar.block(0,0,computed_rank,N2);
				Ac = Mat(N1,computed_rank);
				Ar = Mat(computed_rank,N2);
				return;
			}

			col_bases.push_back(col_index);
			v=col/col(int(row_index));
			row_bases.push_back(int(row_index));
			row = K->getRow(row_indices[row_index], col_indices);
			u	=	row;
			Uvec.push_back(u);
			Vvec.push_back(v);
			// Ac.col(computed_rank) = col;
			// Ar.row(computed_rank) = row;
			AcVec.push_back(col);
			ArVec.push_back(row);
			computed_rank = 1;

			remaining_row_ind.erase(row_index);
			remaining_col_ind.erase(col_index);

			double normS	=	0.0;
			double prev_normS = DBL_MAX;
			this->maxAbsVector(row, remaining_col_ind, max, col_index);

			// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && (abs(prev_normS-normS) >= tol_ACA*prev_normS || abs(prev_normS-normS) >= tol_ACA) && computed_rank < min) {
			while (abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
			// while(u.norm()*v.norm() > tol_ACA*normS &&  computed_rank < min) {
			// while (abs(row(int(col_index))) > tol_ACA && abs(col(int(row_index))) > tol_ACA && abs(prev_normS-normS) >= tol_ACA*prev_normS && computed_rank < min) {
				col_bases.push_back(int(col_index));
			  col = K->getCol(col_indices[col_index], row_indices);
				col_temp = col;
			  for (int l=0; l<=computed_rank-1; ++l){
			    for (size_t s = 0; s < Vvec[l].size(); s++) {
			      col(s)	-=	Uvec[l](col_index)*Vvec[l](s);
			    }
			  }
			  this->maxAbsVector(col, remaining_row_ind, max, row_index);
			  row_bases.push_back(int(row_index));
			  for (size_t l = 0; l < col.size(); l++) {
			    v(l) = col(l)/col(int(row_index));
			  }
			  row = K->getRow(row_indices[row_index], col_indices);
				row_temp = row;
			  for (int l=0; l<=computed_rank-1; ++l){
			    for (size_t s = 0; s < Uvec[l].size(); s++) {
			      row(s)	-=	Vvec[l](row_index)*Uvec[l](s);
			    }
			  }
			  u	=	row;
			  // if(u.norm()< tol_robust || v.norm()< tol_robust) {
			  //   col_bases.pop_back();
			  //   row_bases.pop_back();
			  //   break;
			  // }
			  Uvec.push_back(u);
			  Vvec.push_back(v);
				// Ac.col(computed_rank) = col_temp;
				// Ar.row(computed_rank) = row_temp;
				AcVec.push_back(col_temp);
				ArVec.push_back(row_temp);

			  ++computed_rank;
				remaining_row_ind.erase(row_index);
			  remaining_col_ind.erase(col_index);
			  if (computed_rank != 2)
			    prev_normS = normS;
			  normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
			  for (int l=0; l<=computed_rank-1; ++l){
			    normS	+=	2*abs(Uvec[l].dot(u) * Vvec[l].dot(v));
			  }
			  normS	=	sqrt(normS);
			  this->maxAbsVector(row, remaining_col_ind, max, col_index);
			}
		}
		// Ac_modified = Ac.block(0,0,N1,computed_rank);
		// Ar_modified = Ar.block(0,0,computed_rank,N2);
		Ac = Mat(N1,computed_rank);
		Ar = Mat(computed_rank,N2);
		for (size_t i = 0; i < computed_rank; i++) {
			Ac.col(i) = AcVec[i];
			Ar.row(i) = ArVec[i];
		}
	}
};
