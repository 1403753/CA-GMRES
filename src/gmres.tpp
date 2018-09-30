template <typename ScalarType>
bool gmres<ScalarType>::is_conj_pair(complex_t a, complex_t b) {
	return (a.real() == b.real() && a.imag() == -b.imag() && a.imag() != 0);
}

template <typename ScalarType>
sparse_status_t gmres<ScalarType>::init_gmres(size_t n, const sparse_matrix_t A, const ScalarType *v, ScalarType **H, ScalarType **Q, std::vector<pair_t, mkl_allocator<pair_t>> &theta_vals, size_t s, size_t m) {
	sparse_status_t 			                      stat = SPARSE_STATUS_SUCCESS;
	ScalarType 						                      beta;
	ScalarType                                  *w;
	ScalarType                                  *H_s;
	ScalarType                                  h_ij;
	ScalarType                                  h_jp1j;
	
	/******************************/
	/*  LAPACKE_dhseqr variables  */
	/******************************/
	size_t                                      ilo = 1;        // A hasn't been balanced by ?gebal. Therefore, ilo = 1
	size_t                                      ihi = s;        // A hasn't been balanced by ?gebal. Therefore, ilo = 1
	char                                        job = 'E';      // eigenvalues only are required
	const char                                  compz = 'N';    // no Schur vectors are computed
	ScalarType                                  *wr, *wi;       // real and imag part of ritz values
	ScalarType																	*scale;
	size_t                                      i, j;           // index in for-loops
	
	w = (ScalarType *)mkl_calloc(n, sizeof(ScalarType), 64);if(w == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	H_s = (ScalarType *)mkl_calloc((s + 1)*s, sizeof(ScalarType), 64);if(H_s == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	
	beta = cblas_dnrm2 (n, v, 1);
	
	// Q is in col-major!
	for(i = 0; i < n; ++i)
		(*Q)[i] = v[i] / beta;
	
	if (n < 15) {
		for(i = 0; i < n; ++i) {
			for(size_t k = 0; k < m+1; ++k) {
				std::cout << (*Q)[k*n + i] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	
	// std::cout << std::endl;
	// for(size_t o = 0; o < n; ++o)
		// std::cout << w[o] << "\n";	
	
	for(j = 0; j < s; ++j) {
		mv(A, &(*Q)[j*n], &w, 1);
		for(i = 0; i < j+1; ++i) {
			
			h_ij = cblas_ddot(n, w, 1, &(*Q)[i*n], 1);
			H_s[s*i + j] = h_ij;
			(*H)[m*i + j] = h_ij;
			
			cblas_daxpy (n, h_ij*(-1), &(*Q)[i*n], 1, w, 1);
		}
		h_jp1j = cblas_dnrm2(n, w, 1);
		H_s[s*(j+1) + j] = h_jp1j;
		(*H)[m*(j+1) + j] = h_jp1j;
		
		cblas_daxpy(n, 1 / h_jp1j, w, 1, &(*Q)[(j+1)*n], 1);	
	}

	printf("\n============= H_s:\n");
	for(size_t i = 0; i < s+1; ++i) {
		for(size_t j = 0; j < s; ++j) {
			printf("%2.2f ", H_s[i*s + j]);
		}
		printf("\n");
	}
	
	
	wr = (double *)mkl_malloc(s*sizeof(double), 64);if(wr == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	wi = (double *)mkl_malloc(s*sizeof(double), 64);if(wi == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	scale = (double *)mkl_malloc(s*sizeof(double), 64);if(scale == NULL){return SPARSE_STATUS_ALLOC_FAILED;}



// balancing prob. not needed!
	// job = 'N';
//LAPACKE_dgebal( int matrix_layout, char job, lapack_int n, double* a, lapack_int lda, lapack_int* ilo, lapack_int* ihi, double* scale );
	// LAPACKE_dgebal(LAPACK_ROW_MAJOR, job, s, H_s, s, &ilo, &ihi, scale);	
	
	job = 'E';
	
//LAPACKE_dhseqr(int matrix_layout, char job, char compz, lapack_int n, lapack_int ilo, lapack_int ihi, double *h, lapack_int ldh, double *wr, double *wi, double *z, lapack_int ldz);
	LAPACKE_dhseqr(LAPACK_ROW_MAJOR, job, compz, s, ilo, ihi, H_s, s, wr, wi, nullptr, s);	
	
	std::cout << std::endl;
	for(size_t o = 0; o < s; ++o)
		printf("%.10f,\t%.10f\n", wr[o], wi[o]);	

	modified_leya_ordering(s, wr, wi, theta_vals);
	
	mkl_free(w);
	mkl_free(wr);
	mkl_free(wi);
	mkl_free(H_s);
	mkl_free_buffers();
	
	return stat;
}

template <typename ScalarType>
sparse_status_t gmres<ScalarType>::mv(const sparse_matrix_t A, const ScalarType *x, ScalarType **y, size_t s) {
	sparse_status_t 			stat = SPARSE_STATUS_SUCCESS;
	struct matrix_descr 	descr;
	
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	
	mkl_sparse_set_mv_hint(A, SPARSE_OPERATION_NON_TRANSPOSE, descr, s);
	mkl_sparse_optimize(A);
		
//sparse_status_t mkl_sparse_d_mv (sparse_operation_t operation, double alpha, const sparse_matrix_t A, struct matrix_descr descr, const double *x, double beta, double *y);
	stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descr, x, 0, *y);
	if(stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("MatMult failed");
	
	mkl_free_buffers();
	return stat;
}

template <typename ScalarType>
sparse_status_t gmres<ScalarType>::modified_leya_ordering(size_t s, ScalarType *wr, ScalarType *wi, std::vector<pair_t, mkl_allocator<pair_t>> &theta_vals) {
	sparse_status_t 			stat = SPARSE_STATUS_SUCCESS;
	std::vector<pair_t, mkl_allocator<pair_t>> ritz_vals;


//	ritz_vals.reserve(sizeof(wr) / sizeof(wr[0]));
	ritz_vals.reserve(s);

	for(size_t i = 0; i < s; ++i) {
		ritz_vals.push_back(pair_t(1, complex_t(wr[i], wi[i])));
	}
		
	std::stable_sort(ritz_vals.begin( ), ritz_vals.end( ), [ ]( const pair_t& lhs, const pair_t& rhs ) {
			return lhs.second.real() < rhs.second.real();
	});
	
	for(size_t i = 0; i < ritz_vals.size() - 1;) {
		if(ritz_vals.at(i).second.imag() > 0) {
			if(i + 3 < ritz_vals.size() && ritz_vals.at(i).second == ritz_vals.at(i + 2).second && ritz_vals.at(i + 1).second == ritz_vals.at(i + 3).second) {
				ritz_vals.at(i).first++;
				ritz_vals.erase(ritz_vals.begin() + i + 2);
				ritz_vals.at(i + 1).first++;
				ritz_vals.erase(ritz_vals.begin() + i + 2);
			} else {
				i += 2;
			}
		} else {
			if(ritz_vals.at(i).second == ritz_vals.at(i + 1).second) {
				ritz_vals.at(i).first++;
				ritz_vals.erase(ritz_vals.begin() + i + 1);
			} else {
				i++;
			}
		}
	}

	ritz_vals.shrink_to_fit();
	
	std::vector<size_t, mkl_allocator<size_t>> k_index;

	theta_vals.reserve(ritz_vals.size());
	k_index.reserve(ritz_vals.size());
	
	for(size_t i = 0; i < ritz_vals.size(); ++i)
		k_index.push_back(i);

	/*
	*
	* MODIFIED LEJA ORDERING:
	*
	*	input: ritz_vals (unique shifts with multiplicity mu)
	*				 type: std::vector<size_t mu, complex<double> vals> 
	*
	*	output: theta_vals (in modified Leja order with outlist)
	*					type: std::vector<size_t outlist, complex<double> vals>
	*/

	std::cout << "\n\n\t-\tcompute MODIFIED LEJA ORDERING:\n" << std::endl;
	
	complex_t Capacity_old, Capacity = 1;
	size_t L;


	if(std::abs(ritz_vals.back().second) < std::abs(ritz_vals.front().second)) {
		
		if (ritz_vals.front().second.imag() == 0) {

			theta_vals.push_back(pair_t(std::move(k_index.front()), ritz_vals.front().second));
			k_index.erase(k_index.begin());
			L = 0;
		} else {
			
			theta_vals.push_back(pair_t(std::move(k_index.front()), ritz_vals.front().second));
			theta_vals.push_back(pair_t(std::move(k_index.begin()[1]), ritz_vals.begin()[1].second));

			k_index.erase(k_index.begin());
			k_index.erase(k_index.begin());
			L = 1;
		}
	} else {

		if (ritz_vals.back().second.imag() == 0) {

			theta_vals.push_back(pair_t(std::move(k_index.back()), ritz_vals.back().second));
			k_index.pop_back();
			L = 0;
		} else {
			
			theta_vals.push_back(pair_t(std::move(k_index.end()[-2]), ritz_vals.end()[-2].second));
			theta_vals.push_back(pair_t(std::move(k_index.back()), ritz_vals.back().second));
			
			k_index.pop_back();
			k_index.pop_back();				
			L = 1;
		}
	}
	
	
	while(L < ritz_vals.size() - 1) {

		Capacity_old = Capacity;
		Capacity = 1;
		
		for(size_t j = 0; j < L; ++j) {
			Capacity *= std::pow(std::abs(theta_vals.at(L).second - theta_vals.at(j).second), (double) ritz_vals.at(theta_vals.at(j).first).first / (L+1));
		}
				
		for (auto &z: ritz_vals) {
			z.second /= (Capacity / Capacity_old);
		}
		
		for (auto &t: theta_vals)
			t.second /= (Capacity / Capacity_old);
		
		std::vector<complex_t, mkl_allocator<complex_t>> zprod;
		
		for (auto k: k_index) {
			complex_t prod = 1;
			for(size_t j = 0; j < L+1; ++j) {
				prod *= std::pow(std::abs(ritz_vals.at(k).second - ritz_vals.at(theta_vals.at(j).first).second) / Capacity, ritz_vals.at(theta_vals.at(j).first).first);
			}
			zprod.push_back(prod);
		}
		
		std::vector<complex_t, mkl_allocator<complex_t>>::iterator max_zprod;
		max_zprod = std::max_element(zprod.begin( ), zprod.end( ), [ ]( const complex_t& lhs, const complex_t& rhs ) {
			return std::abs(lhs) < std::abs(rhs);
		});
		
		if (*max_zprod == (complex_t)0 ) {
			throw std::invalid_argument("Product to maximize is zero; either there are multiple shifts, or the product underflowed");
		} else if (std::isinf((*max_zprod).real()) ||	std::isnan((*max_zprod).real()) )
			throw std::invalid_argument("Product to maximize is Inf; must have overflowed");
		
		size_t idx = max_zprod - zprod.begin();
		

		
		if(ritz_vals.at(k_index.at(idx)).second.imag() != 0) {
			
			if(ritz_vals.at(k_index.at(idx)).second.imag() < 0) {
				if (!is_conj_pair(ritz_vals.at(k_index.at(idx) - 1).second, ritz_vals.at(k_index.at(idx)).second)) 
					throw std::invalid_argument( "Input out of order");
				
				theta_vals.push_back(pair_t(std::move(k_index.at(idx - 1)), ritz_vals.at(k_index.at(idx) - 1).second));
				theta_vals.push_back(pair_t(std::move(k_index.at(idx)), ritz_vals.at(k_index.at(idx)).second));			

				k_index.erase(k_index.begin() + idx - 1);
				k_index.erase(k_index.begin() + idx - 1);
				
				L += 2;	
			
			
			} else {
				if (!is_conj_pair(ritz_vals.at(k_index.at(idx)).second, ritz_vals.at(k_index.at(idx + 1)).second))
					throw std::invalid_argument( "Input out of order");
				
				theta_vals.push_back(pair_t(std::move(k_index.at(idx)), ritz_vals.at(k_index.at(idx)).second));
				theta_vals.push_back(pair_t(std::move(k_index.at(idx + 1)), ritz_vals.at(k_index.at(idx + 1)).second));			
				
				k_index.erase(k_index.begin() + idx);
				k_index.erase(k_index.begin() + idx);

				L += 2;
			}
		} else {
		
			theta_vals.push_back(pair_t(k_index.at(idx), ritz_vals.at(k_index.at(idx)).second));
			k_index.erase(k_index.begin() + idx);
		
			L++;
		}
	} //end while
	
	std::cout << std:: endl << "ritz_vals :" << std:: endl;
	for (auto &p: ritz_vals) {
		p.second *= Capacity;
		std::cout << p.second << "   \tmult: " << p.first << std::endl;
	}
	
	std::cout << std:: endl << "theta_vals (FINAL):" << std:: endl;
	for (auto &p: theta_vals) {
		p.second *= Capacity;
		std::cout << p.second << "   \toutlist: " << p.first << std:: endl;
	}
	
	return stat;
}
