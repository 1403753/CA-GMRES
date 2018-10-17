/*
 * gmres_ca.hpp
 *
 *  Created on: 17.10.2018
 *      Author: Robert
 */

#ifndef GMRES_CA_HPP_

#include "includes.hpp"
#include "gmres.hpp"
#include <vector>

#define GMRES_CA_HPP_

template <typename ScalarType>
class gmres_ca {
public:
	gmres_ca();

	virtual ~gmres_ca();
};

#endif /* GMRES_HPP_ */