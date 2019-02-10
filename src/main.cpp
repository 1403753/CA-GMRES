
////////////
//  TODO  //
////////////
/*
	clean up
	add init SpMV time measurement for CA-GMRES
*/

#ifndef MAIN_HPP

#include "KSM.hpp"
#include "GMRES_ca.hpp"
#include "GMRES.hpp"
#include "PCILU0_ca.hpp"
#include "PCNone.hpp"
#include "MmtReader.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include "plot.tpp"

#define MAIN_HPP

int main(int argc, char **args) {
	
	sparse_status_t               stat = SPARSE_STATUS_SUCCESS;

	// std::string fname = "watt1";
	// std::string fname = "e05r0000";
	// std::string fname = "sparse9x9complex";
	// std::string fname = "bcsstk18";
	// std::string fname = "bcsstk08";
	// std::string fname = "xenon2";
	// std::string fname = "bmw7st1";
	// std::string fname = "pwtk"; // not working
	// std::string fname = "goodwin";
	// std::string fname = "dwb512";
	// std::string fname = "1138bus";
	// std::string fname = "nasa4704";
	// std::string fname = "minitest"; // too small, to work properly
	// std::string fname = "CA-ILU(0)"; // too small, to work properly

	if (argc > 5) {
		std::string type = std::string(args[2]);
		
		std::istringstream iss_threads(args[1]);
		size_t numThreads;
		iss_threads >> numThreads;
		mkl_set_num_threads(numThreads);	
		
		if (type == "residuals" && argc == 10) {

			std::string fname = args[3];
			std::istringstream iss_st(args[4]);
			std::istringstream iss_s1(args[5]);
			std::istringstream iss_t1(args[6]);
			std::istringstream iss_s2(args[7]);
			std::istringstream iss_t2(args[8]);

			size_t st;
			size_t s1;
			size_t t1;
			size_t s2;
			size_t t2;
			iss_st >> st;
			iss_s1 >> s1;
			iss_t1 >> t1;
			iss_s2 >> s2;
			iss_t2 >> t2;
			
			std::string title = args[9]; // "\\shortstack{\\footnotesize\\,1e+04\\,$\\times$\\,1e+04\\,diagonal,\\,logspace\\,eigs,\\,cond\\,1e+05\\\\\\footnotesize\\,Residual\\,2-norm,\\,log\\,scale}";
			generate_residual_plot(fname, title, st, s1, t1, s2, t2);

		} else if (type == "speedup" && argc > 7) {
			
			std::istringstream iss_its(args[3]);
			std::istringstream iss_s(args[4]);
			std::istringstream iss_t(args[5]);
			std::string title = std::string(args[6]);
			
			size_t its;
			size_t s;
			size_t t;
			iss_its >> its;
			iss_s >> s;
			iss_t >> t;

			std::vector<std::string> fnames;

			fnames.reserve(argc - 7);

			for (size_t i = 7; i < (size_t) argc; ++i) {
				fnames.push_back(std::string(args[i]));
			}

			generate_speedup_plot(&fnames, title, s, t, its);

		} else if (type == "threads" && argc > 7) {
			
			std::istringstream iss_its(args[3]);
			std::istringstream iss_s(args[4]);
			std::istringstream iss_t(args[5]);
			std::string title = std::string(args[6]);
			
			size_t its;
			size_t s;
			size_t t;
			iss_its >> its;
			iss_s >> s;
			iss_t >> t;
			
			std::string fname = std::string(args[7]);

			generate_thread_plot(fname, title, s, t, its);

		}

	}

	return stat;
}

#endif