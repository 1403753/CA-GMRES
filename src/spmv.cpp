/*
 * spmv.cpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */
#include "spmv.hpp"
#include <iostream>
#include <exception>

typedef double ScalarType;

spmv::spmv() {
	
}

//void spmv::mv(std::vector< std::map<unsigned int, double>> A, const std::vector<double> x, std::vector<double> res) {
//void spmv::mv(const std::vector<ScalarType> A, const std::vector<ScalarType> x, std::vector<ScalarType> &res) {
void spmv::mv(const sparse_matrix_t A, const double *x, double **y) {
/*	
	{
		std::vector<std::map<unsigned int, ScalarType>> _A_(5);
		
		_A_[0][0] = 1;
		_A_[1][1] = 1;
		_A_[2][2] = 1;
		_A_[3][3] = 1;
		_A_[4][4] = 1;
		
		viennacl::vector<ScalarType> vcl_res(res->size());
		viennacl::vector<ScalarType> vcl_x(x.size());

	//  viennacl::coordinate_matrix<ScalarType> vcl_A(res.size(), x.size());
		viennacl::compressed_matrix<ScalarType> vcl_A(5, 5);
	//  viennacl::matrix<ScalarType> vcl_A(res.size(), x.size());
		
		// initialize viennacl data
		copy(x.begin(), x.end(), vcl_x.begin());
		copy(_A_, vcl_A);
		
		
		for (size_t i = 0; i < vcl_A.size1(); ++i) {
			for (size_t j = 0; j < vcl_A.size2(); ++j)
			std::cout << vcl_A(i, j) << ' ';
		std::cout << std::endl;
		}
		

		vcl_res = viennacl::linalg::prod(vcl_A, vcl_x);
		
		copy(vcl_res.begin(), vcl_res.end(), res->begin());
		
		std::cout << "res: ";
		for(auto s : vcl_res)
			std::cout << s << ' ';
		
		std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
	}
	*/
	
	///////MKL
	{		
//		const size_t m = 5;
//		const size_t n = 3;
//		double x[m] = {1, 1, 1, 1, 1};
/*
		MKL_INT rows_ptr[] = {0, 3, 5, 8, 11, 13};
		MKL_INT col_indx[] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4};
		double values[] = 	 {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
*/
		sparse_status_t stat;		

/*
		MKL_INT rows_ptr[] = {0, 1, 2, 3, 4, 5};
		MKL_INT col_indx[] = {0, 0, 1, 0, 0};
		double values[] = 	 {1, 2, 1, 1, 1};
		
		sparse_matrix_t A;
		
//	initialize A
//	stat = mkl_sparse_d_create_csr (A, indexing, rows, cols, rows_start, rows_end, col_indx, values)
    stat = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, m, n, rows_ptr, rows_ptr + 1, col_indx, values);
	*/	
		
		struct matrix_descr descr;
		descr.type = SPARSE_MATRIX_TYPE_GENERAL;
		
		//sparse_status_t mkl_sparse_d_mv (sparse_operation_t operation, double alpha, const sparse_matrix_t A, struct matrix_descr descr, const double *x, double beta, double *y);
		stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descr, &x[0], 0, *y);
		

		
		
		if(stat != SPARSE_STATUS_SUCCESS)
			throw std::invalid_argument("MatMult failed");
		
		mkl_free_buffers();
	}
	
}


spmv::~spmv() {
	
}

