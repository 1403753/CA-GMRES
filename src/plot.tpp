
sparse_status_t res_solve(std::string fname, std::string fdir, std::string title, size_t n, sparse_matrix_t *A_mkl, IKSPType *kspType, IPCType *pcType, double *b, size_t s, size_t t, size_t st) {
	
	std::ofstream file;
	sparse_status_t               stat;
	
	double                 		   	rTol = 1e-11;                   // the relative (possibly preconditioned) residual norm || A*x_{k + 1} - b || / || A*x_0 - b ||
                                                                // == || r_{k+1} || / || r_0 ||
	double                    		aTol = 1e-50;                   // the absolute (possibly preconditioned) residual norm || A*x_{k + 1} - b || == || r_{k+1} ||
	double                        dTol = 1e+4;                    // the divergence tolerance, amount (possibly preconditioned) residual norm can increase
	size_t                        maxit = 1000;                   // maximum number of iterations to use	
	double *x;
	KSP_ ksp;	
	
	x = (double *) mkl_calloc(n, sizeof(double), 64);if(x == NULL){return SPARSE_STATUS_ALLOC_FAILED;}

	ksp.setOptions(rTol, aTol, dTol, maxit, true);

	ksp.setOperator(A_mkl);

	ksp.setKSPType(kspType);
	ksp.setPCType(pcType);

	ksp.setUp();

	stat = ksp.solve(x, b);

	std::vector<std::pair<size_t, double>>* rHist = ksp.getRHist();

	file.open(fdir);

	if(dynamic_cast<GMRES*>(kspType)) {
		file << fname << " " << st << " " << title << std::endl;
	} else if (dynamic_cast<GMRES_ca*>(kspType)) {
		std::streamsize ss = std::cout.precision();
		file << fname << " " << s << " " << t << " " << std::setprecision(2) << std::scientific << ((GMRES_ca*) kspType)->getRcondMin() << " " << ((GMRES_ca*) kspType)->getRcondMax() << std::setprecision(ss) << std::endl;
	} else exit(1);
	
	for (auto h:*rHist) {
		file << h.first << " " << std::scientific << h.second << std::endl;
	}

	mkl_free(x);

	file.close();
	
	return stat;
}

sparse_status_t generate_residual_plot(std::string fname, std::string title, size_t st, size_t s1, size_t t1, size_t s2, size_t t2) {

	sparse_status_t               stat;
	sparse_matrix_t               A_mkl;                          // n x n matrix A
	double                        *b;
	double                        *tx;
	double                        *r;
	PCILU0_ca                     ilu0;                           // PCType
	PCNone                        pcnone;                         // PCType
	MmtReader											mmtReader;
	GMRES_ca											gmres_ca(s1, t1, NEWTON);     // KSPType
	GMRES                         gmres(st);                      // KSPType
	
	sparse_index_base_t           indexing;	
	size_t                        n, m;
	size_t                        *rows_start;
	size_t                        *rows_end;
	size_t                        *col_indx;
	double                        *values;
	struct matrix_descr           descr;
	std::ofstream                 file;
	std::ifstream                 ifile;

	// read matrix
	stat = mmtReader.read_matrix_from_file(std::string("../matrix_market/") + fname + std::string(".mtx"), &A_mkl);

	descr.type = SPARSE_MATRIX_TYPE_GENERAL;

	stat = mkl_sparse_d_export_csr(A_mkl, &indexing, &n, &m, &rows_start, &rows_end, &col_indx, &values); // n is needed

	b = (double *) mkl_malloc(n * sizeof(double), 64);if(b == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	tx = (double *) mkl_malloc(n * sizeof(double), 64);if(tx == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	r = (double *) mkl_malloc(n * sizeof(double), 64);if(r == NULL){return SPARSE_STATUS_ALLOC_FAILED;}

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, time(NULL)); // Seed with time	

		
	// file.open("tx.vec_bcsskt18");
	
	// for (size_t i = 0; i < n; ++i) {
		// tx[i] = (gsl_ran_flat(rng, -1, 1) + std::sin(2*M_PI*i/n)) *2;
		// // tx[i] = 1;		
		// file << std::scientific << std::setprecision(20) << tx[i] << std::endl;
	// }
	
	// file.close();


	// start with initial guess
	// for (size_t i = 0; i < n; ++i)
		// x[i] = 0.1*(i%3 + 1);
	
	ifile.open("tx.vec");
	for (size_t i = 0; i < n; ++i) {
		ifile >> tx[i];
	}
	ifile.close();

	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A_mkl, descr, tx, 0, b);

	res_solve(fname, "../gnuplot/gmres_newt_small_s.dat", title, n, &A_mkl, &gmres_ca, &pcnone, b, s1, t1, s1*t1);
	gmres_ca.setBasis(MONOMIAL);
	res_solve(fname, "../gnuplot/gmres_mono_small_s.dat", title, n, &A_mkl, &gmres_ca, &pcnone, b, s1, t1, s1*t1);

	gmres_ca.setS(s2);
	gmres_ca.setT(t2);

	res_solve(fname, "../gnuplot/gmres_mono_large_s.dat", title, n, &A_mkl, &gmres_ca, &pcnone, b, s2, t2, s2*t2);	
	gmres_ca.setBasis(NEWTON);
	res_solve(fname, "../gnuplot/gmres_newt_large_s.dat", title, n, &A_mkl, &gmres_ca, &pcnone, b, s2, t2, s2*t2);
	res_solve(fname, "../gnuplot/gmres_standard.dat", title, n, &A_mkl, &gmres, &pcnone, b, st, st, st);

	gsl_rng_free(rng);

	stat = mkl_sparse_destroy(A_mkl);
	mkl_free(tx);
	mkl_free(b);
	mkl_free(r);
	mkl_free_buffers();

	return stat;
}

sparse_status_t speedup_solve(std::string fname, std::string title, size_t idx, KSP_ *ksp, size_t s, size_t t, size_t its) {
	
	sparse_matrix_t               A_mkl;                          // n x n matrix A
	sparse_index_base_t           indexing;
	GMRES_ca											gmres_ca(s, t, NEWTON);         // KSPType
	GMRES                         gmres(s*t);                     // KSPType	
	sparse_status_t               stat;	
	size_t                        n, m;
	size_t                        *rows_start;
	size_t                        *rows_end;
	size_t                        *col_indx;
	double                        *values;
	double                        *b;
	double                        *x;
	double                        *tx;	
	std::ofstream                 file;
	struct matrix_descr           descr;
	float                         gmres_ca_rtime;
	float                         avg_gmres_ca_rtime = 0;
	float                         avg_SDO = 0;
	float                         avg_BCGS = 0;
	float                         avg_TSQR = 0;
	float                         avg_SpMV_ca = 0;
	float                         avg_INIT = 0;
	float                         avg_MGS = 0;
	float                         avg_SpMV = 0;
	
	MmtReader											mmtReader;

	stat = mmtReader.read_matrix_from_file(std::string("../matrix_market/") + fname + std::string(".mtx"), &A_mkl);
		
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	
	stat = mkl_sparse_d_export_csr(A_mkl, &indexing, &n, &m, &rows_start, &rows_end, &col_indx, &values); // n is needed

	b = (double *) mkl_malloc(n * sizeof(double), 64);if(b == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	x = (double *) mkl_malloc(n * sizeof(double), 64);if(x == NULL){return SPARSE_STATUS_ALLOC_FAILED;}	
	tx = (double *) mkl_malloc(n * sizeof(double), 64);if(tx == NULL){return SPARSE_STATUS_ALLOC_FAILED;}

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, time(NULL)); // Seed with time	

	for (size_t i = 0; i < n; ++i) {
		tx[i] = gsl_ran_flat(rng, -1, 1) + std::sin(2*M_PI*i/n);
		// tx[i] = 1;		
	}	
	
	stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A_mkl, descr, tx, 0, b);	

	stat = ksp->setOperator(&A_mkl);

	for (size_t i = 0; i < its; ++i) {
		
		std::cout << "iter: " << i + 1 << std::endl;
		
		stat = ksp->setKSPType(&gmres_ca);
		stat = ksp->setUp();

		std::fill(x, x + n, 0);		
		stat = ksp->solve(x, b);	

		stat = ksp->setKSPType(&gmres);

		std::fill(x, x + n, 0);
		stat = ksp->solve(x, b);

		gmres_ca_rtime = gmres_ca.getSDO() + gmres_ca.getBCGS() + gmres_ca.getTSQR() + gmres_ca.getSpMV() + gmres_ca.getINIT();
		avg_gmres_ca_rtime += gmres_ca_rtime/its;

		avg_SDO += gmres_ca.getSDO()/its;
		avg_BCGS += gmres_ca.getBCGS()/its;
		avg_TSQR += gmres_ca.getTSQR()/its;
		avg_SpMV_ca += gmres_ca.getSpMV()/its;
		avg_INIT += gmres_ca.getINIT()/its;

		avg_MGS += gmres.getMGS()/its;
		avg_SpMV += gmres.getSpMV()/its;

	}

		file.open("../gnuplot/gmres.dat", std::ios_base::app);
		file << idx << " " << std::string("\\\\footnotesize\\\\,") + fname << " " << avg_MGS / avg_gmres_ca_rtime << " " << avg_SpMV / avg_gmres_ca_rtime << std::endl;
		file.close();

		file.open("../gnuplot/gmres_ca.dat", std::ios_base::app);
		file << idx << " " << avg_SDO / avg_gmres_ca_rtime << " " << avg_BCGS / avg_gmres_ca_rtime << " " << avg_TSQR / avg_gmres_ca_rtime \
			<< " " << avg_SpMV_ca / avg_gmres_ca_rtime << " " << avg_INIT / avg_gmres_ca_rtime << std::endl;
		file.close();

	mkl_sparse_destroy(A_mkl);
	gsl_rng_free(rng);
	mkl_free(x);
	mkl_free(tx);
	mkl_free(b);

	return stat;
}

sparse_status_t generate_speedup_plot(std::vector<std::string> *fnames, std::string title, size_t s, size_t t, size_t its){
	
	double                 		   	rTol = 1e-11;                   // the relative (possibly preconditioned) residual norm || A*x_{k + 1} - b || / || A*x_0 - b ||
																																// == || r_{k+1} || / || r_0 ||
	double                    		aTol = 1e-50;                   // the absolute (possibly preconditioned) residual norm || A*x_{k + 1} - b || == || r_{k+1} ||
	double                        dTol = 1e+4;                    // the divergence tolerance, amount (possibly preconditioned) residual norm can increase
	size_t                        maxit = 1000;                   // maximum number of iterations to use
	sparse_status_t               stat;
	KSP_                          ksp;							 						  // linear solver context	
	PCILU0_ca                     ilu0;                           // PCType
	PCNone                        pcnone;                         // PCType
	std::ofstream                 file;	

	// file.open("../gnuplot/gmres_ca.dat");
		// file << "nr" << " " << "SDO" << " " << "BCGS" << " " << "SpMV" << std::endl;
	// file.close();
	
	// file.open("../gnuplot/gmres.dat");
		// file << "nr" << " " << "title" << " " << "MGS" << " " << "SpMV" << " " << title << std::endl;
	// file.close();

	stat = ksp.setOptions(rTol, aTol, dTol, maxit, true);
	stat = ksp.setPCType(&pcnone);
	
	for (size_t i = 0; i < fnames->size(); ++i) {
		stat = speedup_solve(fnames->at(i), title, i+1, &ksp, s, t, its);
	}

	mkl_free_buffers();

	return stat;
}