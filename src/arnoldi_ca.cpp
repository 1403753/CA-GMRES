/*
 * arnoldi_ca.cpp
 *
 *  Created on: 11.05.2018
 *      Author: Robert
 */

#include "arnoldi_ca.hpp"

arnoldi_ca::arnoldi_ca() {
	// Constructor
	
}

bool is_conj_pair(std::complex<double> a, std::complex<double> b) {
	return (a.real() == b.real() && a.imag() == -b.imag() && a.imag() != 0);
}

void arnoldi_ca::print() {
	
	float rtime, ptime, mflops;
	long long flpops;
	float irtime, iptime, imflops;
  long long iflpops;
	
	std::cout.precision(3);
	
	/*
	*
	*	compute ritz values used in arnoldi(s)
	*
	*/
	
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, time(NULL));

	const size_t n = 14;
	
	
	
	size_t alpha = 1;
	const char job = 'E';
	const char compz = 'N';
	double *H, *wr, *wi;
	wr = (double *)mkl_malloc(n*sizeof(double), 64);
	wi = (double *)mkl_malloc(n*sizeof(double), 64);

	H = (double *)mkl_calloc(n*n, sizeof(double), 64);

	for(size_t i = 0; i < n; ++i) 
		for(size_t j = i; j < n; ++j) {
			H[i*n + j] = gsl_rng_uniform_int(rng, 100) + gsl_rng_uniform(rng);
			if(i + 1 < n)
				H[(i + 1)*n + j] = gsl_rng_uniform_int(rng, 100) + gsl_rng_uniform(rng);
		}
	/*
	for(size_t i = 0; i < n; ++i) {
		for(size_t j = 0; j < n; ++j)
			std::cout << H[i*n + j] << ' ';
		std::cout << std::endl;
	}
	*/
	for(size_t i = 0; i < n; ++i) {
		for(size_t j = 0; j < n; ++j)
			std::cout << H[i*n + j] << ' ';
		std::cout << std::endl;
	}

	/*
	*	compute ritz values of H
	*/
//LAPACKE_dhseqr(int matrix_layout, char job, char compz, lapack_int n, lapack_int ilo, lapack_int ihi, double *h, lapack_int ldh, double *wr, double *wi, double *z, lapack_int ldz);
	LAPACKE_dhseqr(LAPACK_ROW_MAJOR, job, compz, n, alpha, n, H, n, wr, wi, nullptr, n);
	
	std::cout << "H (after):" << std::endl;

	for(size_t i = 0; i < n; ++i) {
		for(size_t j = 0; j < n; ++j)
			std::cout << H[i*n + j] << ' ';
		std::cout << std::endl;
	}
	
	
	std::vector<std::pair<size_t, std::complex<double>>> ritz_vals;


//	ritz_vals.reserve(sizeof(wr) / sizeof(wr[0]));
	ritz_vals.reserve(n);

	for(size_t i = 0; i < n; ++i) {
		ritz_vals.push_back(std::pair<size_t, std::complex<double>>(1, std::complex<double>(wr[i], wi[i])));
	}
		
	std::stable_sort(ritz_vals.begin( ), ritz_vals.end( ), [ ]( const std::pair<size_t, std::complex<double>>& lhs, const std::pair<size_t, std::complex<double>>& rhs ) {
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
	
	std::vector<std::pair<size_t, std::complex<double>>> theta_vals;
	std::vector<size_t> k_index;

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
	
	std::complex<double> C_old, C = 1;
	size_t L;


	if(std::abs(ritz_vals.back().second) < std::abs(ritz_vals.front().second)) {
		
		if (ritz_vals.front().second.imag() == 0) {

			theta_vals.push_back(std::pair<size_t, std::complex<double>>(std::move(k_index.front()), ritz_vals.front().second));
			k_index.erase(k_index.begin());
			L = 0;
		} else {
			
			theta_vals.push_back(std::pair<size_t, std::complex<double>>(std::move(k_index.front()), ritz_vals.front().second));
			theta_vals.push_back(std::pair<size_t, std::complex<double>>(std::move(k_index.begin()[1]), ritz_vals.begin()[1].second));

			k_index.erase(k_index.begin());
			k_index.erase(k_index.begin());
			L = 1;
		}
	} else {

		if (ritz_vals.back().second.imag() == 0) {

			theta_vals.push_back(std::pair<size_t, std::complex<double>>(std::move(k_index.back()), ritz_vals.back().second));
			k_index.pop_back();
			L = 0;
		} else {
			
			theta_vals.push_back(std::pair<size_t, std::complex<double>>(std::move(k_index.end()[-2]), ritz_vals.end()[-2].second));
			theta_vals.push_back(std::pair<size_t, std::complex<double>>(std::move(k_index.back()), ritz_vals.back().second));
			
			k_index.pop_back();
			k_index.pop_back();				
			L = 1;
		}
	}
	
	
	while(L < ritz_vals.size() - 1) {

		C_old = C;
		C = 1;
		
		for(size_t j = 0; j < L; ++j) {
			C *= std::pow(std::abs(theta_vals.at(L).second - theta_vals.at(j).second), (double) ritz_vals.at(theta_vals.at(j).first).first / (L+1));
		}
				
		for (auto &z: ritz_vals) {
			z.second /= (C / C_old);
		}
		
		for (auto &t: theta_vals)
			t.second /= (C / C_old);
		
		std::vector<std::complex<double>> zprod;
		
		for (auto k: k_index) {
			std::complex<double> prod = 1;
			for(size_t j = 0; j < L+1; ++j) {
				prod *= std::pow(std::abs(ritz_vals.at(k).second - ritz_vals.at(theta_vals.at(j).first).second) / C, ritz_vals.at(theta_vals.at(j).first).first);
			}
			zprod.push_back(prod);
		}
		
		std::vector<std::complex<double>>::iterator max_zprod;
		max_zprod = std::max_element(zprod.begin( ), zprod.end( ), [ ]( const std::complex<double>& lhs, const std::complex<double>& rhs ) {
			return std::abs(lhs) < std::abs(rhs);
		});
		
		if (*max_zprod == (std::complex<double>)0 ) {
			throw std::invalid_argument("Product to maximize is zero; either there are multiple shifts, or the product underflowed");
		} else if (std::isinf((*max_zprod).real()) ||	std::isnan((*max_zprod).real()) )
			throw std::invalid_argument("Product to maximize is Inf; must have overflowed");
		
		size_t idx = max_zprod - zprod.begin();
		

		
		if(ritz_vals.at(k_index.at(idx)).second.imag() != 0) {
			
			if(ritz_vals.at(k_index.at(idx)).second.imag() < 0) {
				if (!is_conj_pair(ritz_vals.at(k_index.at(idx) - 1).second, ritz_vals.at(k_index.at(idx)).second)) 
					throw std::invalid_argument( "Input out of order");
				
				theta_vals.push_back(std::pair<size_t, std::complex<double>>(std::move(k_index.at(idx - 1)), ritz_vals.at(k_index.at(idx) - 1).second));
				theta_vals.push_back(std::pair<size_t, std::complex<double>>(std::move(k_index.at(idx)), ritz_vals.at(k_index.at(idx)).second));			

				k_index.erase(k_index.begin() + idx - 1);
				k_index.erase(k_index.begin() + idx - 1);
				
				L += 2;	
			
			
			} else {
				if (!is_conj_pair(ritz_vals.at(k_index.at(idx)).second, ritz_vals.at(k_index.at(idx + 1)).second))
					throw std::invalid_argument( "Input out of order");
				
				theta_vals.push_back(std::pair<size_t, std::complex<double>>(std::move(k_index.at(idx)), ritz_vals.at(k_index.at(idx)).second));
				theta_vals.push_back(std::pair<size_t, std::complex<double>>(std::move(k_index.at(idx + 1)), ritz_vals.at(k_index.at(idx + 1)).second));			
				
				k_index.erase(k_index.begin() + idx);
				k_index.erase(k_index.begin() + idx);

				L += 2;
			}
		} else {
		
			theta_vals.push_back(std::pair<size_t, std::complex<double>>(k_index.at(idx), ritz_vals.at(k_index.at(idx)).second));
			k_index.erase(k_index.begin() + idx);
		
			L++;
		}
	} //end while
	
	if (n < 15) {
		std::cout << std:: endl << "ritz_vals :" << std:: endl;
		for (auto p: ritz_vals) {
			std::cout << p.second*C << "\tmult: " << p.first << std::endl;
		}
		
		std::cout << std:: endl << "theta_vals (FINAL):" << std:: endl;
		for (auto p: theta_vals) {
			std::cout << p.second*C << "\toutlist: " << p.first << std:: endl;
		}
	}

	
	
	
	/*
	*
	* Apply MODIFIED GIVENS ROTATION to H (used in CA-GMRES)
	*
	*/
	std::cout << "\n\n==============================GIVENS" << std::endl;

	
	
	
	// re-initialize H	
	for(size_t i = 0; i < n; ++i) 
		for(size_t j = 0; j < n; ++j) {
			if(j < i-1)
				H[i*n + j] = 0;
			else
				H[i*n + j] = gsl_rng_uniform_int(rng, 100) + gsl_rng_uniform(rng);
		}
	
	if(n < 2) {
		std::cout << std::endl << "H:\n";

		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j)
				std::cout << H[n*i + j] << "\t";
			std::cout << std::endl;
		}
	}	
	
	double c = 0;
	double s = 0;
	double a = 0;
	double b = 0;
	
	if (PAPI_flops(&irtime, &iptime, &iflpops, &imflops) != PAPI_OK)
		exit(1);
	
	for (size_t i = 0; i < n - 1; ++i) {
			a = H[i*n + i];
			b = H[(i + 1)*n + i];
// void cblas_drotg (double *a, double *b, double *c, double *s);
			cblas_drotg(&a, &b, &c, &s);

// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
			cblas_drot(n - i, &H[i*n + i], 1, &H[(i + 1)*n + i], 1, c, s);
	}
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) != PAPI_OK)
		exit(1);
	
	PAPI_shutdown();
	std::cout << n << ", " << rtime << std::endl;	
	
	if(n < 2) {
		std::cout << std::endl << std::endl << "H:\n";

		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j)
				std::cout << H[n*i + j] << "\t";
			std::cout << std::endl;
		}
	}
	
	
	
	
	// re-initialize H	
	for(size_t i = 0; i < n; ++i) 
		for(size_t j = 0; j < n; ++j) {
			if(j < i-1)
				H[i*n + j] = 0;
			else
				H[i*n + j] = gsl_rng_uniform_int(rng, 100) + gsl_rng_uniform(rng);
		}
	
	if(n < 2) {
		std::cout << std::endl << "H:\n";

		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j)
				std::cout << H[n*i + j] << "\t";
			std::cout << std::endl;
		}
	}
	
	double d1 = INT32_MAX;
	double d2 = INT32_MAX;
	double x1 = 0;
	double y1 = 0;
	double param[5]{};
	
	if (PAPI_flops(&irtime, &iptime, &iflpops, &imflops) != PAPI_OK)
			exit(1);
		
	for (size_t i = 0; i < n - 1; ++i) {
		//std::cout << std::endl << "d1: " << d1 << ", d2: "<< d2;
		x1 = H[i*n + i];
		y1 = H[(i + 1)*n + i];
//	void cblas_drotmg (double *d1, double *d2, double *x1, const double y1, double *param);
		cblas_drotmg (&d1, &d2, &x1, y1, param);
	
//	void cblas_drotm (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double *param);
		cblas_drotm (n - i, &H[i*n + i], 1, &H[(i + 1)*n + i], 1, param);

	}
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) != PAPI_OK)
		exit(1);
	
	PAPI_shutdown();
	std::cout << n << ", " << rtime << std::endl;	
	
	if(n < 2) {
		std::cout << std::endl << std::endl << "H:\n";

		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j)
				std::cout << H[n*i + j] << "\t";
			std::cout << std::endl;
		}
	}
	
	std::cout << std::endl << "d1: " << d1 << ", d2: "<< d2 <<  std::endl << std::endl;
	std::cout << "param: ";
	for (size_t j = 0; j < 5; ++j)
		std::cout << param[j] << ' ';
	
	std::cout << std::endl << std::endl;	
	
	
	
	
	
	
	mkl_free(H);
	mkl_free(wr);
	mkl_free(wi);
	gsl_rng_free(rng);
	mkl_free_buffers();
}

arnoldi_ca::~arnoldi_ca() {
	// Destructor
	
}



	/* NOT NEEDED!!
	HAD TO ITERATE OVER RITZ VALUES AND SWAP COMPLEX NUMBERS BECAUSE SORTING ALGO
	WAS UNSTABLE. REPLACEMENT WITH STABLE VERSION MADE THIS CODE OBSOLETE.
	
	for (size_t i = 0; i < n; ++i)
		if(std::imag(ritz_vals.at(i)) != 0) {
			if(std::imag(ritz_vals.at(i)) < 0) {
				std::swap(ritz_vals[i], ritz_vals[i + 1]);
				++i;
			} else ++i;
			
			if (std::real(ritz_vals.at(i-1)) != std::real(ritz_vals.at(i))) {
				std::cout << "error!\n";	
				exit(-1);
			}
		}
	std::cout << "ritz:\n";
	for(auto r: ritz_vals)
		std::cout << r << "\n";
	
	
	NOT NEEDED!!
	NORMAL GIVENS ROTATION (Valgrind error) AND MODIFIED GIVENS ROTATION EXAMPLES WITH 3X3 MATRICES
	
	{
		const int m = 3;
		double vec[m*m] = {9,0,2,
											 0,7,3,
											 2,3,1};
		double c = 0;
		double s = 0;
		double a = vec[0];
		double b = vec[m*2];

		for (size_t j = 0; j < m - 1; ++j) {
			for (size_t i = m; --i > j;) {
				a = vec[j*m + j];
				b = vec[i*m + j];
// void cblas_drotg (double *a, double *b, double *c, double *s);
				cblas_drotg(&a, &b, &c, &s);

// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
				cblas_drot(m - j, &vec[j*m + j], 1, &vec[i*m + j], 1, c, s);
			}
		}
		
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < m; ++j)
				std::cout << vec[m*i + j] << ' ';
			std::cout << std::endl;
		}
	
		std::cout << std::endl << "c: " << c << ", s: "<< s <<  std::endl << std::endl;
		double vec2[m*m] = {9,3,2,
												2,7,3,
												0,3,1};
		double d1 = 1;
		double d2 = 1;
		double x1;
		double param[5]{};
		
		for (size_t i = 0; i < m - 1; ++i) {
			x1 = vec2[i*m + i];
//	void cblas_drotmg (double *d1, double *d2, double *x1, const double y1, double *param);
			cblas_drotmg (&d1, &d2, &x1, vec2[(i + 1)*m + i], param);
	
//	void cblas_drotm (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double *param);
			cblas_drotm (m - i, &vec2[i*m + i], 1, &vec2[(i + 1)*m + i], 1, param);
		}
	
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < m; ++j)
				std::cout << vec2[m*i + j] << ' ';
			std::cout << std::endl;
		}
	
		std::cout << std::endl << "d1: " << d1 << ", d2: "<< d2 <<  std::endl << std::endl;
		std::cout << "param: ";
		for (size_t j = 0; j < 5; ++j)
			std::cout << param[j] << ' ';
		
		std::cout << std::endl;
	}
	*/
	
	/*
	*
	* TEST DATA for modified Leja order
	*
	*/
	/*
	double wr[n] = {0.3795668652,  0.3795668652, -0.3795668652, -0.3795668652, -1, 0.3414450976, -0.3795668652, -0.3795668652, -1, .01, 0};
	double wi[n] = {10.7694800092, -10.7694800092,  0.7694800092, -0.7694800092,  0, 					 	0,  0.7694800092, -0.7694800092,  0,   0, 0};

	double wr[n] = {2, 2, 9, 8};
	double wi[n] = {1, -1, 0, 0};

	double wr[n] = {-0.23278561593838409394,
								-0.23278561593838409394,
								-0.17965204298588813292,
								-0.17965204298588813292,
								0.00000000000000000000,
								0.00000000000000000000,
								0.00000000000000004290,
								1.00000000000000155431,
								1.46557123187676840992,
								2.00000000000000000000,
								2.35930408597177665442};
	double wi[n] = {0.79255199251544805605,
								 - 0.79255199251544805605,
									0.90301314585700431792,
									- 0.90301314585700431792,
									0.00000000000000000000,
									0.00000000000000000000,
									0.00000000000000000000,
									0.00000000000000000000,
									0.00000000000000000000,
									0.00000000000000000000,
									0.00000000000000000000};

	double wr[n] = {-18.26253971751562232839,
									-18.26253971751562232839,
									-0.76880864595760334268,
									7.33095544939625387570,
									7.33095544939625387570,
									23.03908177957729108698,
									49.47271159472406765190,
									75.45757786431740044009,
									110.84124849021880265809,
									127.88031559350544341669,
									182.32163992436568378253};

	double wi[n] = {8.56750903964367971355,
									-8.56750903964367971355,
									0.00000000000000000000,
									46.23141328511955805425,
									-46.23141328511955805425,
									0,0,0,0,0,0};
	*/