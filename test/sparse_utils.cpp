#include "sparse_utils.hpp"

// taken from igraph ../src/structural_properties.c
sparse_status_t neighborhood(const Mtx_CSR *A,  // CSR MATRIX
												std::vector<size_t *> &v_In, // input set
												std::vector<size_t *> &v_Res, // result vertices for some level of order
												size_t order, // the reachabilty order
												Gr_Part gPart
								)
{
  
  size_t n = A->n; // calculate number of nodes in the graph
	size_t *Ap = A->row_ptr;
	size_t *Ai = A->col_indx;
	
  std::deque<size_t> q; // double-ended queue (C++: std::deque)
	//* igraph_dqueue_pop :=  Remove the head. (C++: pop_front)
	//* igraph_dqueue_pop_back :=  Remove the tail. (C++: pop_back)
	//* igraph_dqueue_push := Append an element to the end of the queue (head side) (C++: push_back) (there is no push_back)

	//* igraph_dqueue_back := The last element in the queue (C++: back)
	//* igraph_dqueue_head := The first element in the queue (C++: front)
	//* igraph_dqueue_size := number of elements in the queue (C++: size)
	
  size_t *added; // added vertices
  std::vector<size_t> neis; // neighbors of a vertex

  if (order < 0) {
    return SPARSE_STATUS_INVALID_VALUE;
  }
  
	// allocate memory for added nodes with size of "number of all nodes in the graph" of type "int"
	added = (size_t*) mkl_calloc(n, sizeof(size_t), 64); if(added == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}

	for(std::vector<size_t*>::reverse_iterator it = v_In.rbegin(); it != v_In.rend(); ++it) { //iterate over input set. (for each vid do .... )
		ptrdiff_t node = std::distance(Ap, *it); // get the current vertex_id, aka node_id;
    added[node]=**it+1; // add node to added vertices-vector in location "node (vertex_id)" and set it to vertexnumber + 1
		v_Res.push_back(Ap + node); // add current vertex to result.
    
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
					if (added[nei] != **it+1) { // only add if the added vertices-vector in location "nei (vertex_id)" is not the vertexnumber + 1
																	 //* in other words this should mean if the neighbor is not the current vertex in the input set or has not been added yet!
																	 //* *added is kind of like the marker for "has been visited"
						added[nei]=**it + 1; // mark as visited
						q.push_back(nei); // add the vertex to the dqueue
						q.push_back(actdist + 1); // add its distance to the queue
						v_Res.push_back(Ap + nei); // add the vertex to the result 
					}
				}
      } else {
				/* we just count them but don't add them to q */ // just add them to the result, but not process them in the dqueue
				for (size_t j=0; j<nneis; j++) { // iterate over all neighbors
					size_t nei= neis.at(j); // get current neighbor
					if (added[nei] != **it+1) {
						added[nei] = **it + 1;
						v_Res.push_back(Ap + nei); // add the vertex to the result 
					}
				}
      }
    } /* while q not empty */
		
  } //* end for every vertex in the input set

	if (v_In.size() > 1) {
		std::sort( v_Res.begin(), v_Res.end() );
		v_Res.erase( std::unique( v_Res.begin(), v_Res.end() ), v_Res.end() );
	}

	mkl_free(added);
  return SPARSE_STATUS_SUCCESS;
}






sparse_status_t my_permute(const Mtx_CSR *A, Mtx_CSR *dest, const size_t n, const size_t *pinv, const size_t *q) {

	size_t t, j, k, inz = 0, *Ap, *Ai, *Cp, *Ci;
	double *Cx, *Ax;

	Ap = A->row_ptr;
	Ai = A->col_indx;
	Ax = A->values;
	
	Cp = dest->row_ptr;
	Ci = dest->col_indx;
	Cx = dest->values;

	for (k = 0 ; k < n ; k++) {
			Cp[k] = inz;                   /* row k of C is row pinv[k] of A */
			j = pinv ? (pinv[k]) : k ;
			for (t = Ap[j] ; t < Ap[j+1] ; t++) {
					if (Cx) Cx[inz] = Ax[t] ;  /* column i of A is column q[i] of C */
					Ci[inz++] = q ? (q[Ai [t]]) : Ai[t] ;
			}
	}
	Cp[n] = inz;                       /* finalize the last row of C */
			
	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t at_plus_a(
	const size_t n,      /* number of columns in matrix A. */
	const size_t nz,     /* number of nonzeros in matrix A */
	size_t *colptr,      /* column pointer of size n+1 for matrix A. */
	size_t *rowind,      /* row indices of size nz for matrix A. */
	size_t *bnz,         /* out - on exit, returns the actual number of
														 nonzeros in matrix A'*A. */
	size_t **b_colptr,   /* out - size n+1 */
	size_t **b_rowind    /* out - size *bnz */
	  )
{
/*
 * Purpose
 * =======
 *
 * Form the structure of A'+A. A is an n-by-n matrix in column oriented
 * format represented by (colptr, rowind). The output A'+A is in column
 * oriented format (symmetrically, also row oriented), represented by
 * (b_colptr, b_rowind).
 *
 */
	size_t i, j, k, col, num_nz;
	size_t *t_colptr, *t_rowind; /* a column oriented form of T = A' */
	size_t *marker;

	marker = (size_t*) mkl_malloc(n * sizeof(size_t), 64); if(marker == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	t_colptr = (size_t*) mkl_malloc((n+1) * sizeof(size_t), 64); if(t_colptr == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	t_rowind = (size_t*) mkl_malloc(nz * sizeof(size_t), 64); if(t_rowind == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	
	/* Get counts of each column of T, and set up column pointers */
	for (i = 0; i < n; ++i) marker[i] = 0;
	for (j = 0; j < n; ++j) {
		for (i = colptr[j]; i < colptr[j+1]; ++i)
			++marker[rowind[i]];
	}
	t_colptr[0] = 0;
	for (i = 0; i < n; ++i) {
		t_colptr[i+1] = t_colptr[i] + marker[i];
		marker[i] = t_colptr[i];
	}

	/* Transpose the matrix from A to T */
	for (j = 0; j < n; ++j)
		for (i = colptr[j]; i < colptr[j+1]; ++i) {
			col = rowind[i];
			t_rowind[marker[col]] = j;
			++marker[col];
		}


	/* ----------------------------------------------------------------
		 compute B = A + T, where column j of B is:

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

		/* Add pattern of column A_*k to B_*j */
		for (i = colptr[j]; i < colptr[j+1]; ++i) {
			k = rowind[i];
			if ( marker[k] != j ) {
				marker[k] = j;
				++num_nz;
			}
		}

		/* Add pattern of column T_*k to B_*j */
		for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
			k = t_rowind[i];
			if ( marker[k] != j ) {
				marker[k] = j;
				++num_nz;
			}
		}
	}
		
	*bnz = num_nz;
  
	/* Allocate storage for A+A' */
	*b_colptr = (size_t*) mkl_malloc((n+1) * sizeof(size_t), 64); if(b_colptr == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	if ( *bnz) {
		*b_rowind = (size_t*) mkl_malloc(*bnz * sizeof(size_t), 64); if(b_rowind == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	}
    
	/* Zero the diagonal flag */
	for (i = 0; i < n; ++i) marker[i] = -1;
    
	/* Compute each column of B, one at a time */
	num_nz = 0;
	for (j = 0; j < n; ++j) {
		(*b_colptr)[j] = num_nz;
	
		/* Flag the diagonal so it's not included in the B matrix */
		marker[j] = j;

		/* Add pattern of column A_*k to B_*j */
		for (i = colptr[j]; i < colptr[j+1]; ++i) {
				k = rowind[i];
				if ( marker[k] != j ) {
			marker[k] = j;
			(*b_rowind)[num_nz++] = k;
				}
		}

		/* Add pattern of column T_*k to B_*j */
		for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
			k = t_rowind[i];
			if ( marker[k] != j ) {
				marker[k] = j;
				(*b_rowind)[num_nz++] = k;
			}
		}
	}
	(*b_colptr)[n] = num_nz;
		 
	mkl_free(marker);
	mkl_free(t_colptr);
	mkl_free(t_rowind);	
	
	return SPARSE_STATUS_SUCCESS;
}