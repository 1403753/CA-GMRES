template <typename ScalarType>
sparse_status_t gmres<ScalarType>::tsqr(ScalarType *V, ScalarType *Q, ScalarType *R_, const size_t m, const size_t n, const size_t st) {

	size_t lwork = -1;
	size_t lmwork = -1;
	size_t tsize = -1;
	size_t info = 0;
	size_t i, j;
	ScalarType *work, *t, *tquery, *workquery, *mwork, *mworkquery;
	char side = 'L', trans = 'N';
	sparse_status_t stat = SPARSE_STATUS_SUCCESS;
		
	for(i = 0; i < n; ++i) 
		Q[i*m + i] = 1;

	tquery = (ScalarType *)mkl_malloc(5 * sizeof(ScalarType), 64);if (tquery == NULL)throw std::invalid_argument("malloc error: allocating memory for tquery failed.");
	workquery = (ScalarType *)mkl_malloc(sizeof(ScalarType), 64);if (workquery == NULL)throw std::invalid_argument("malloc error: allocating memory for workquery failed.");
	
	/* Workspace query */
	dgeqr(&m, &n, V, &m, tquery, &tsize, workquery, &lwork, &info);
	
	tsize = tquery[0];
	mkl_free(tquery);
	
	t = (ScalarType *)mkl_malloc(tsize * sizeof(ScalarType), 64);if (t == NULL)throw std::invalid_argument("malloc error: allocating memory for t failed.");
	
	lwork = workquery[0];
	mkl_free(workquery);	
	
	work = (ScalarType *)mkl_malloc(lwork * sizeof(ScalarType), 64);if (work == NULL)throw std::invalid_argument("malloc error: allocating memory for work failed.");
	
	dgeqr(&m, &n, V, &m, t, &tsize, work, &lwork, &info);
	mkl_free(work);
	
	mworkquery = (ScalarType *)mkl_malloc(sizeof(ScalarType), 64);if (mworkquery == NULL)throw std::invalid_argument("malloc error: allocating memory for mworkquery failed.");
	
	/************************/
	/*  store Q explicitly  */
	/************************/
	
	/* Workspace query */
	dgemqr(&side, &trans, &m, &n, &n, V, &m, t, &tsize, Q, &m, mworkquery, &lmwork, &info);

	lmwork = mworkquery[0];
	mkl_free(mworkquery);
		
	mwork = (ScalarType *)mkl_malloc(lmwork * sizeof(ScalarType), 64);if (mwork == NULL)throw std::invalid_argument("malloc error: allocating memory for mwork failed.");
	
	dgemqr(&side, &trans, &m, &n, &n, V, &m, t, &tsize, Q, &m, mwork, &lmwork, &info);
	mkl_free(mwork);

	mkl_free(t);
	
	for(i = 0; i < n; ++i) {
		for(j = i; j < n; ++j) {
			R_[(st + n + 1)*j + i] = V[m*j + i];
		}
	}

	if (info < 0) stat = SPARSE_STATUS_EXECUTION_FAILED;

	return stat;
}

	// std::cout << "AAAAARRR:" << std::endl;
	// for(i = 0; i < n; ++i) {
		// for(j = i; j < n; ++j) {
			// std::cout << R_[(st + n + 1)*j + i] << " ";
		// }
		// std::cout << std::endl;
	// }
	// printf("\n============= R:\n");
	// for(size_t i = 0; i < n; ++i) {
		// for(size_t j = 0; j < n; ++j) {
			// printf("%f, ", V[n*j + i]);
		// }
		// printf("\n");
	// }
	// if (m < 15) {
		// printf("\n============= Q:\n");
		// for(size_t i = 0; i < m; ++i) {
			// for(size_t j = 0; j < n; ++j) {
				// printf("%f, ", Q[m*j + i]);
			// }
			// printf("\n");
		// }
	// }
	// std::cout << "\n\n";	
	// for(size_t j = 0; j < n*m; ++j)
		// printf("%f, ", V[j]);
	// std::cout << "\n";

	// if (m < 15) {
		// printf("\n============= V:\n");
		// for(size_t i = 0; i < m; ++i) {
			// for(size_t j = 0; j < n; ++j) {
				// std::cout << V[m*j + i] << " ";
			// }
			// std::cout << std::endl;
		// }
	// }	
	

	// cblas_dtrmm (CblasColMajor, CblasRight, CblasUpper, CblasNoTrans,  CblasNonUnit, 
								// m, n, 1, V, m, Q, m);
	// if (m < 15) {
		// printf("\n============= QR:\n");
		// for(size_t i = 0; i < m; ++i) {
			// for(size_t j = 0; j < n; ++j) {
				// printf("%.20f, ", Q[m*j + i]);
			// }
			// printf("\n");
		// }
	// }
	// ScalarType *C;
	// C = (ScalarType *)mkl_calloc(n*n, sizeof(ScalarType), 64);if (C == NULL)throw std::invalid_argument("malloc error: allocating memory for C failed.");
	
	// cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, n, n, m, 1, Q, m, Q, m, 0, C, n);
	
	// printf("\n============= C:\n");
	// for(size_t i = 0; i < n; ++i) {
		// for(size_t j = 0; j < n; ++j) {
			// printf("%.20f, ", C[n*j + i]);
		// }
		// printf("\n");
	// }
	
	// mkl_free(C);