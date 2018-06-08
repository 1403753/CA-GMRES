/*
 * spmv.hpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#ifndef SPMV_HPP_
#include <vector>
#include <map>

#include "mkl_spblas.h"
#include "mkl.h"

/*
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/prod.hpp"  
*/
#define SPMV_HPP_

class spmv {
public:
	spmv();
	static void mv(const std::vector<double>, const std::vector<double>, std::vector<double>&);
	virtual ~spmv();
};

#endif /* SPMV_HPP_ */
