/*
 * gmres.cpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */
#include "gmres.hpp"

template <typename ScalarType>
gmres<ScalarType>::gmres() {
	
}

template <typename ScalarType>
gmres<ScalarType>::~gmres() {
	
}

/*	NOT NEEDED

	VIENNACL SpMV:
	{

	#include "viennacl/vector.hpp"
	#include "viennacl/compressed_matrix.hpp"
	#include "viennacl/matrix.hpp"
	#include "viennacl/linalg/prod.hpp"  

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