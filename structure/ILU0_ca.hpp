#ifndef ILU0_CA_HPP
#define ILU0_CA_HPP

#include "IPCType.hpp"

// what happens:
	// compute structure at_plus_a(A_mtx) ->output: b_rowptr and partition A with metis 
	// sort A_mtx with amml(s)
	// permute A_mtx -> col_indx out of order
	// create A_mkl and sort col_indx (maybe create own algorithm)
	// export A_mtx from A_mkl again
	// factorize A_mtx to get ILU(0)-values
	// find subsets alpha_p p == #threads depending on matrix size
	// openmp: find s-step dependencies for each alpha_p -> beta_p -> gamma_p -> delta_p 

class ILU0_ca : public IPCType{
public:
	ILU0_ca();
	sparse_status_t find(Mtx_CSR *A_mtx,
											 std::vector<size_t>& alphaArr, 
											 std::vector<std::vector<size_t>>& betaArr,
											 std::vector<std::vector<size_t>>& gammaArr,
											 std::vector<std::vector<size_t>>& deltaArr);
	sparse_status_t setUp();
};

#endif


	////////////////////
	//  METIS OUTPUT  //
	////////////////////

	// size_t k = -1;
	// std::cout << "\nMetis K-Way:" << std::endl;
	// while (++k != (size_t) mkl_get_max_threads()) {
		// for(size_t i = 0; i < n; i++) {
			// if (part[i] == k) {
				// std::cout << i << ": " << part[i] << std::endl;
				// break;
			// }
		// }
	// }


	/////////////////////
	//  PAPI TEMPLATE  //
	/////////////////////

	// float                         rtime, ptime, mflops;
	// long long                     flpops;		
	// if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		// exit(1);
	// // -> function here <-
	// if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		// exit(1);
	// PAPI_shutdown();
	// std::cout << "runtime: " << rtime << " seconds." << std::endl;