#include "sparse_utils.hpp"

sparse_status_t extract(const Mtx_CSR *A,  // CSR input-Matrix. Column indices have to be in increasing order!
												Mtx_CSR *dest, // extracted CSR output-Matrix
												std::vector<size_t> &rows, // input rows
												std::vector<size_t> &cols) // input columns
{
	
	if(!std::is_sorted(cols.begin(), cols.end()))
		std::sort(cols.begin(), cols.end());

	size_t n = A->n; // number of nodes in the original graph
	size_t *Ap = A->row_ptr;
	size_t *Ai = A->col_indx;
	double *Ax = A->values;
	
	size_t *t_row_ptr;
	size_t *t_col_indx;
	double *t_values;
		
	size_t nrows = rows.size();
	
	t_row_ptr = (size_t*) mkl_malloc((nrows + 1) * sizeof(size_t), 64); if(t_row_ptr == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	t_col_indx = (size_t*) mkl_malloc(nrows*n * sizeof(size_t), 64); if(t_col_indx == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	t_values = (double*) mkl_malloc(nrows*n * sizeof(double), 64); if(t_values == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}

	t_row_ptr[0] = 0;

	size_t nz = 0;

	for (size_t i = 0; i < nrows; ++i) {
		size_t optimizer = 0;
		for (size_t j = Ap[rows[i]]; j < Ap[rows[i]+1]; ++j) {
			auto result = std::find(cols.begin() + optimizer, cols.end(), Ai[j]);
			if(result != cols.end()) {				
				t_col_indx[nz] = result - cols.begin();
				t_values[nz++] = Ax[j];
				++optimizer;
			}
			t_row_ptr[i+1] = nz;
		}
	}
	
	mkl_realloc(t_col_indx, nz);
	mkl_realloc(t_values, nz);
		
	dest->n = nrows;
	dest->row_ptr = t_row_ptr;
	dest->col_indx = t_col_indx;
	dest->values = t_values;
	
	return SPARSE_STATUS_SUCCESS;
}

//! taken from igraph ../src/structural_properties.c 
sparse_status_t neighborhood(const Mtx_CSR *A,  // CSR Matrix
												std::vector<size_t> &v_In, // input set
												std::vector<size_t> &v_Res, // result vertices for some level of order
												size_t order, // the reachabilty order
												Gr_Part gPart) // the part of the graph to operate on
{
  
  size_t n = A->n; // number of nodes in the graph
	size_t *Ap = A->row_ptr;
	size_t *Ai = A->col_indx;
	v_Res.clear();
	v_Res.reserve(100);
  std::deque<size_t> q; // double-ended queue (C++: std::deque)
	//* igraph_dqueue_pop :=  Remove the head. (C++: pop_front)
	//* igraph_dqueue_pop_back :=  Remove the tail. (C++: pop_back)
	//* igraph_dqueue_push := Append an element to the end of the queue (head side) (C++: push_back) (there is no push_back)

	//* igraph_dqueue_back := The last element in the queue (C++: back)
	//* igraph_dqueue_head := The first element in the queue (C++: front)
	//* igraph_dqueue_size := number of elements in the queue (C++: size)
	
  size_t *marker; // added vertices
  std::vector<size_t> neis; // neighbors of a vertex

  if (order < 0) {
    return SPARSE_STATUS_INVALID_VALUE;
  }
  
	// allocate memory for added nodes with size of "number of all nodes in the graph" of type "int"
	marker = (size_t*) mkl_calloc(n, sizeof(size_t), 64); if(marker == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}

	for(std::vector<size_t>::reverse_iterator it = v_In.rbegin(); it != v_In.rend(); ++it) { //iterate over input set. (for each vid do .... )
		size_t node = *it; // get the current vertex_id, aka node_id;
    marker[node]=*it + 1; // add node to added vertices-vector in location "node (vertex_id)" and set it to vertexnumber + 1
		v_Res.push_back(node); // add current vertex to result.
    
		if (order > 0) {
			q.push_back(node); // add vertex to the end of the dqueue if order > 0
			q.push_back(0); // add 0 afterwards to the end of the dqueue (this number is the distance of the node to itself)
    }

    while(!q.empty()) { // as long as dqueue is not empty do ...
			size_t actnode = q.front(); // remove the head (the node)
			q.pop_front();
			size_t actdist = q.front(); // remove the head (the distance/order of the node)
			q.pop_front();
      size_t nneis; // declare nneis := number of neighbors
			for(size_t i = *(Ap + actnode); i < *(Ap + actnode + 1); ++i) {
				switch(gPart) {
					case GRAPH_LOWER:
						if(actnode > Ai[i])
							neis.push_back(Ai[i]); // get all neighbors (outdegree) from the current vertex and store in neis
						break;
					case GRAPH_UPPER:
						if(actnode <= Ai[i])
							neis.push_back(Ai[i]); // get all neighbors (outdegree) from the current vertex and store in neis
						break;
					case GRAPH_COMPLETE:
						neis.push_back(Ai[i]); // get all neighbors (outdegree) from the current vertex and store in neis
						break;
					default:
						return SPARSE_STATUS_INVALID_VALUE;
				}
			}
			
			nneis = neis.size(); // calculate size of neighbors
      
      if (actdist<order-1) { // ask if the node-distance is smaller than the order (in first iteration actdist == 0 and order-1 == 0 for order==1)
				/* we add them to the q */
				for (size_t j = 0; j < nneis; j++) { // iterate over all neighbors
					size_t nei= neis.at(j); // get current neighbor
					if (marker[nei] != *it+1) { // only add if the added vertices-vector in location "nei (vertex_id)" is not the vertexnumber + 1
																	 //* in other words this should mean if the neighbor is not the current vertex in the input set or has not been added yet!
																	 //* *marker is kind of like "has been visited"
						marker[nei]=*it + 1; // mark as visited
						q.push_back(nei); // add the vertex to the dqueue
						q.push_back(actdist + 1); // add its distance to the queue
						v_Res.push_back(nei); // add the vertex to the result 
					}
				}
      } else {
				/* we just count them but don't add them to q */ // just add them to the result, but not process them in the dqueue
				for (size_t j=0; j<nneis; j++) { // iterate over all neighbors
					size_t nei= neis.at(j); // get current neighbor
					if (marker[nei] != *it+1) {
						marker[nei] = *it + 1;
						v_Res.push_back(nei); // add the vertex to the result 
					}
				}
      }
    } /* while q not empty */
		
  } //* end for every vertex in the input set

	if (v_In.size() > 1) {
		std::sort(v_Res.begin(), v_Res.end());
		v_Res.erase(std::unique( v_Res.begin(), v_Res.end() ), v_Res.end());
	}

	mkl_free(marker);
  return SPARSE_STATUS_SUCCESS;
}



/*! \file taken from SuperLU_MT_3.1/SRC/get_perm_c.c and modified
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy)
*/
sparse_status_t permute_Mtx(const Mtx_CSR *A, Mtx_CSR *dest, const size_t *pinv, const size_t *q) {

	size_t t, j, k, n = A->n, inz = 0, *Ap, *Ai, *Cp, *Ci;
	double *Cx, *Ax;

	Ap = A->row_ptr;
	Ai = A->col_indx;
	Ax = A->values;
	
	Cp = dest->row_ptr;
	Ci = dest->col_indx;
	Cx = dest->values;

	for (k = 0; k < n; k++) {
			Cp[k] = inz; // row k of C is row pinv[k] of A
			j = pinv ? (pinv[k]) : k;
			for (t = Ap[j]; t < Ap[j+1]; t++) {
					if (Cx)
						Cx[inz] = Ax[t]; // column i of A is column q[i] of C
					Ci[inz++] = q ? (q[Ai[t]]) : Ai[t];
			}
	}
	Cp[n] = inz; // finalize the last row of C
			
	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t at_plus_a(const Mtx_CSR *A,
	size_t *bnz,         /* out - on exit, returns the actual number of nonzeros in matrix A'*A. */
	size_t **b_row_ptr,  /* out - size n+1 */
	size_t **b_col_indx  /* out - size *bnz */
) {
/*
 * Purpose
 * =======
 *
 * Form the structure of A'+A. A is an n-by-n matrix in row oriented
 * format represented by (Ap, Ai). The output A'+A is in row
 * oriented format (symmetrically, also column oriented), represented by
 * (b_row_ptr, b_col_indx).
 *
 */
	size_t i, j, k, row, num_nz;
	size_t n = A->n;
	size_t *Ap = A->row_ptr;
	size_t *Ai = A->col_indx;
	size_t nz = Ap[n];

	size_t *t_row_ptr, *t_col_indx; /* a row oriented form of T = A' */
	size_t *marker;

	marker = (size_t*) mkl_malloc(n * sizeof(size_t), 64); if(marker == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	t_row_ptr = (size_t*) mkl_malloc((n+1) * sizeof(size_t), 64); if(t_row_ptr == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	t_col_indx = (size_t*) mkl_malloc(nz * sizeof(size_t), 64); if(t_col_indx == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	
	/* Get counts of each row of T, and set up row pointers */
	for (i = 0; i < n; ++i) marker[i] = 0;
	for (j = 0; j < n; ++j) {
		for (i = Ap[j]; i < Ap[j+1]; ++i)
			++marker[Ai[i]];
	}
	t_row_ptr[0] = 0;
	for (i = 0; i < n; ++i) {
		t_row_ptr[i+1] = t_row_ptr[i] + marker[i];
		marker[i] = t_row_ptr[i];
	}

	/* Transpose the matrix from A to T */
	for (j = 0; j < n; ++j)
		for (i = Ap[j]; i < Ap[j+1]; ++i) {
			row = Ai[i];
			t_col_indx[marker[row]] = j;
			++marker[row];
		}


	/* ----------------------------------------------------------------
		 compute B = A + T, where row j of B is:

		 Struct (B_*j) = Struct (A_*k) UNION Struct (T_*k)

		 do not include the diagonal entry
		 ---------------------------------------------------------------- */

	/* Zero the diagonal flag */
	for (i = 0; i < n; ++i) marker[i] = -1;

	/* First pass determines number of nonzeros in B */
	num_nz = 0;
	for (j = 0; j < n; ++j) {
		/* Flag the diagonal so it's not included in the B matrix */
		marker[j] = j;

		/* Add pattern of row A_*k to B_*j */
		for (i = Ap[j]; i < Ap[j+1]; ++i) {
			k = Ai[i];
			if ( marker[k] != j ) {
				marker[k] = j;
				++num_nz;
			}
		}

		/* Add pattern of row T_*k to B_*j */
		for (i = t_row_ptr[j]; i < t_row_ptr[j+1]; ++i) {
			k = t_col_indx[i];
			if ( marker[k] != j ) {
				marker[k] = j;
				++num_nz;
			}
		}
	}
		
	*bnz = num_nz;
  
	/* Allocate storage for A+A' */
	*b_row_ptr = (size_t*) mkl_malloc((n+1) * sizeof(size_t), 64); if(b_row_ptr == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	if ( *bnz) {
		*b_col_indx = (size_t*) mkl_malloc(*bnz * sizeof(size_t), 64); if(b_col_indx == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	}
    
	/* Zero the diagonal flag */
	for (i = 0; i < n; ++i) marker[i] = -1;
    
	/* Compute each row of B, one at a time */
	num_nz = 0;
	for (j = 0; j < n; ++j) {
		(*b_row_ptr)[j] = num_nz;
	
		/* Flag the diagonal so it's not included in the B matrix */
		marker[j] = j;

		/* Add pattern of row A_*k to B_*j */
		for (i = Ap[j]; i < Ap[j+1]; ++i) {
				k = Ai[i];
				if ( marker[k] != j ) {
			marker[k] = j;
			(*b_col_indx)[num_nz++] = k;
				}
		}

		/* Add pattern of row T_*k to B_*j */
		for (i = t_row_ptr[j]; i < t_row_ptr[j+1]; ++i) {
			k = t_col_indx[i];
			if ( marker[k] != j ) {
				marker[k] = j;
				(*b_col_indx)[num_nz++] = k;
			}
		}
	}
	(*b_row_ptr)[n] = num_nz;
		 
	mkl_free(marker);
	mkl_free(t_row_ptr);
	mkl_free(t_col_indx);	
	
	return SPARSE_STATUS_SUCCESS;
}