sparse_status_t GMRES_ca::tsqr(double *V, double *Q, double *R_, const size_t m, const size_t n, const size_t st) {

	size_t lwork = -1;
	size_t lmwork = -1;
	size_t tsize = -1;
	size_t info = 0;
	size_t i, j;
	double *work, *t, *tquery, *workquery, *mwork, *mworkquery;
	char side = 'L', trans = 'N';
	sparse_status_t stat = SPARSE_STATUS_SUCCESS;
		
	for(i = 0; i < n; ++i) 
		Q[i*m + i] = 1;

	tquery = (double *)mkl_malloc(5 * sizeof(double), 64);if (tquery == NULL)throw std::invalid_argument("malloc error: allocating memory for tquery failed.");
	workquery = (double *)mkl_malloc(sizeof(double), 64);if (workquery == NULL)throw std::invalid_argument("malloc error: allocating memory for workquery failed.");
	
	/* Workspace query */
	dgeqr(&m, &n, V, &m, tquery, &tsize, workquery, &lwork, &info);
	
	tsize = tquery[0];
	mkl_free(tquery);
	
	t = (double *)mkl_malloc(tsize * sizeof(double), 64);if (t == NULL)throw std::invalid_argument("malloc error: allocating memory for t failed.");
	
	lwork = workquery[0];
	mkl_free(workquery);	
	
	work = (double *)mkl_malloc(lwork * sizeof(double), 64);if (work == NULL)throw std::invalid_argument("malloc error: allocating memory for work failed.");
	
	dgeqr(&m, &n, V, &m, t, &tsize, work, &lwork, &info);
	mkl_free(work);
	
	mworkquery = (double *)mkl_malloc(sizeof(double), 64);if (mworkquery == NULL)throw std::invalid_argument("malloc error: allocating memory for mworkquery failed.");
	
	/************************/
	/*  store Q explicitly  */
	/************************/
	
	/* Workspace query */
	dgemqr(&side, &trans, &m, &n, &n, V, &m, t, &tsize, Q, &m, mworkquery, &lmwork, &info);

	lmwork = mworkquery[0];
	mkl_free(mworkquery);
		
	mwork = (double *)mkl_malloc(lmwork * sizeof(double), 64);if (mwork == NULL)throw std::invalid_argument("malloc error: allocating memory for mwork failed.");
	
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