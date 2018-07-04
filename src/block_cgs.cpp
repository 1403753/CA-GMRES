/*
 * block_cgs.cpp
 *
 *  Created on: 07.05.2018
 *      Author: Robert
 *
 */

#include "block_cgs.hpp"

block_cgs::block_cgs() {
	// Constructor
	
}
	
void block_cgs::cgs_tsqr() {
	/*
	* Algorithm 12 (p. 104) implementation of Block CGS
	*	
	*	Input: 
	*	- M blocks of n by m blocked vectors V_k, k = 1,...,M.
	*	
	*	Output: 
	*	- othonormal Q_k with same dimensions as V_k.
	*	- blocked upper triangular R_{ij} (m_i by m_j).
	*/
	
	const size_t M = 10; // total number of V-blocks
	const size_t B = 5;  // total number of vectors in a block
	
	#pragma omp parallel for if(M*B > 10) default(none) \
	schedule(auto) collapse(1)
	for (size_t k = 0; k < M*B; k += B) {
		// TODO: implement compute_qr(k);
	}
}
	
block_cgs::~block_cgs() {
	// Destructor
	
}

	/* NOT NEEDED!!!
		Batch matrix-matrix multiplication for small matrices, not what I'm looking for..
		std::cout << "format: 'MKL_COMPACT_SSE'" << mkl_get_format_compact();
		link: https://software.intel.com/en-us/mkl-developer-reference-c-mkl-get-size-compact
		MKL_INT mkl_dget_size_compact (MKL_INT ld, MKL_INT sd, MKL_COMPACT_PACK format, MKL_INT nm);
	*/
