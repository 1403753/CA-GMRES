template <typename ScalarType>
bool gmres<ScalarType>::is_conj_pair(complex_t a, complex_t b) {
	return (a.real() == b.real() && a.imag() == -b.imag() && a.imag() != 0);
}

template <typename ScalarType>
sparse_status_t gmres<ScalarType>::modified_leya_ordering(size_t s, ScalarType *wr, ScalarType *wi, std::vector<pair_t, mkl_allocator<pair_t>> &theta_vals) {
	sparse_status_t                             stat = SPARSE_STATUS_SUCCESS;
	std::vector<pair_t, mkl_allocator<pair_t>>  ritz_vals;
	size_t                                      i, j, real_idx;

//	ritz_vals.reserve(sizeof(wr) / sizeof(wr[0]));
	ritz_vals.reserve(s);
	
	for(i = 0; i < s; ++i) {
		if(wi[i] == 0) real_idx = i;
	}
	
	for(i = 0; i < s; ++i) {
		if (i != real_idx)
			ritz_vals.push_back(pair_t(1, complex_t(wr[i], wi[i])));		
	}
	
	if (ritz_vals.back().second.imag() > 0.0) {
		ritz_vals.pop_back();
		ritz_vals.push_back(pair_t(1, complex_t(wr[real_idx], wi[real_idx])));		
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
	
	for(i = 0; i < ritz_vals.size(); ++i)
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
		
		for(j = 0; j < L; ++j) {
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
			for(j = 0; j < L+1; ++j) {
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
	
	// std::cout << std:: endl << "ritz_vals :" << std:: endl;
	for (auto &p: ritz_vals) {
		p.second *= Capacity;
		// std::cout << p.second << "   \tmult: " << p.first << std::endl;
	}
	
	// std::cout << std:: endl << "theta_vals (FINAL):" << std:: endl;
	for (auto &p: theta_vals) {
		p.second *= Capacity;
		// std::cout << p.second << "   \toutlist: " << p.first << std:: endl;
	}
	
	return stat;
}
