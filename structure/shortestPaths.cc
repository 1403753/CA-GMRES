/**
 * \ingroup structural
 * \function igraph_shortest_paths
 * \brief The length of the shortest paths between vertices.
 *
 * \param graph The graph object.
 * \param res The result of the calculation, a matrix. A pointer to an
 *        initialized matrix, to be more precise. The matrix will be
 *        resized if needed. It will have the same
 *        number of rows as the length of the \c from
 *        argument, and its number of columns is the number of
 *        vertices in the \c to argument. One row of the matrix shows the
 *        distances from/to a given vertex to the ones in \c to.
 *        For the unreachable vertices IGRAPH_INFINITY is returned.
 * \param from Vector of the vertex ids for which the path length
 *        calculations are done.
 * \param to Vector of the vertex ids to which the path length 
 *        calculations are done. It is not allowed to have duplicated
 *        vertex ids here.
 * \param mode The type of shortest paths to be used for the
 *        calculation in directed graphs. Possible values: 
 *        \clist
 *        \cli IGRAPH_OUT 
 *          the lengths of the outgoing paths are calculated. 
 *        \cli IGRAPH_IN 
 *          the lengths of the incoming paths are calculated. 
 *        \cli IGRAPH_ALL 
 *          the directed graph is considered as an undirected one for
 *          the computation. 
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM 
 *           not enough memory for temporary
 *           data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE 
 *           invalid mode argument.
 *        \endclist
 * 
 * Time complexity: O(n(|V|+|E|)),
 * n is the 
 * number of vertices to calculate, |V| and
 * |E| are the number of vertices and
 * edges in the graph. 
 *
 * \sa \ref igraph_get_shortest_paths() to get the paths themselves, 
 * \ref igraph_shortest_paths_dijkstra() for the weighted version.
 */

int igraph_shortest_paths(const igraph_t *graph,
													igraph_matrix_t *res, 
													const igraph_vs_t from,
													const igraph_vs_t to,
													igraph_neimode_t mode)
{

  long int no_of_nodes=igraph_vcount(graph); //number of nodes of complete graph
  long int no_of_from; //number of nodes of start set?
  long int no_of_to; //number of nodes of result?
  long int *already_counted; // already counted check-array?
  igraph_adjlist_t adjlist; // adjacent list
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL; // double-ended queue (C++: std::deque)
  igraph_vector_int_t *neis; // neighbors of a vertex
  igraph_bool_t all_to; // ???

  long int i, j;
  igraph_vit_t fromvit, tovit; // vertex iterators
  igraph_real_t my_infinity=IGRAPH_INFINITY; //infinity /unreachable value
  igraph_vector_t indexv; // ???

  // if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      // mode != IGRAPH_ALL) {
    // IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
  // }

  // IGRAPH_CHECK(igraph_vit_create(graph, from, &fromvit)); //create vertex iterator
  // IGRAPH_FINALLY(igraph_vit_destroy, &fromvit); //create vertex iterator
  no_of_from=IGRAPH_VIT_SIZE(fromvit); // calculate number of vertices of start set

  // IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode)); // initialize adj. list (mode should be only out!)
  // IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist); // destroy eventually

  already_counted=igraph_Calloc(no_of_nodes, long int); // allocate and initialize with zeros for allready counted vertices
  // if (already_counted==0) { //checks malloc
    // IGRAPH_ERROR("shortest paths failed", IGRAPH_ENOMEM);
  // }
  // IGRAPH_FINALLY(free, already_counted); // destroy memory eventually
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100); // initialize dequeue

  if ( (all_to=igraph_vs_is_all(&to)) ) { // check if all vertices have been initialized in to
    no_of_to=no_of_nodes;
  } else {
    IGRAPH_VECTOR_INIT_FINALLY(&indexv, no_of_nodes);
    IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
    IGRAPH_FINALLY(igraph_vit_destroy, &tovit);
    no_of_to=IGRAPH_VIT_SIZE(tovit);
    for (i=0; !IGRAPH_VIT_END(tovit); IGRAPH_VIT_NEXT(tovit)) {
      long int v=IGRAPH_VIT_GET(tovit);
      if (VECTOR(indexv)[v]) {
	IGRAPH_ERROR("Duplicate vertices in `to', this is not allowed", 
		     IGRAPH_EINVAL);
      }
      VECTOR(indexv)[v] = ++i;
    }
  }

  IGRAPH_CHECK(igraph_matrix_resize(res, no_of_from, no_of_to));
  igraph_matrix_fill(res, my_infinity);

  for (IGRAPH_VIT_RESET(fromvit), i=0; 
       !IGRAPH_VIT_END(fromvit); 
       IGRAPH_VIT_NEXT(fromvit), i++) {
    long int reached=0;
    IGRAPH_CHECK(igraph_dqueue_push(&q, IGRAPH_VIT_GET(fromvit)));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    already_counted[ (long int) IGRAPH_VIT_GET(fromvit) ] = i+1;
    
    IGRAPH_ALLOW_INTERRUPTION();

    while (!igraph_dqueue_empty(&q)) {
      long int act=(long int) igraph_dqueue_pop(&q);
      long int actdist=(long int) igraph_dqueue_pop(&q);

      if (all_to) {
	MATRIX(*res, i, act)=actdist;
      } else {
	if (VECTOR(indexv)[act]) {
	  MATRIX(*res, i, (long int)(VECTOR(indexv)[act]-1)) = actdist;
	  reached++;
	  if (reached==no_of_to) {
	    igraph_dqueue_clear(&q);
	    break;
	  }
	}
      }
      
      neis = igraph_adjlist_get(&adjlist, act);
      for (j=0; j<igraph_vector_int_size(neis); j++) {
        long int neighbor=(long int) VECTOR(*neis)[j];
        if (already_counted[neighbor] == i+1) { continue; }
        already_counted[neighbor] = i+1;
        IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
      }
    }
  }

  /* Clean */
  if (!all_to) {
    igraph_vit_destroy(&tovit);
    igraph_vector_destroy(&indexv);
    IGRAPH_FINALLY_CLEAN(2);
  }

  igraph_Free(already_counted);
  igraph_dqueue_destroy(&q);
  igraph_vit_destroy(&fromvit);
  igraph_adjlist_destroy(&adjlist);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}