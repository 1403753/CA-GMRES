// taken from igraph ../src/structural_properties.c
/** 
 * \function igraph_neighborhood
 * Calculate the neighborhood of vertices.
 * 
 * The neighborhood of a given order of a vertex includes all vertices
 * which are closer to the vertex than the order. Ie. order 0 is
 * always the vertex itself, order 1 is the vertex plus its immediate
 * neighbors, order 2 is order 1 plus the immediate neighbors of the
 * vertices in order 1, etc. 
 * 
 * </para><para> This function calculates the vertices within the
 * neighborhood of the specified vertices.
 * \param graph The input graph.
 * \param res An initialized pointer vector. Note that the objects
 *    (pointers) in the vector will \em not be freed, but the pointer
 *    vector will be resized as needed. The result of the calculation
 *    will be stored here in \c vector_t objects.
 * \param vids The vertices for which the calculation is performed.
 * \param order Integer giving the order of the neighborhood.
 * \param mode Specifies how to use the direction of the edges if a
 *   directed graph is analyzed. For \c IGRAPH_OUT only the outgoing
 *   edges are followed, so all vertices reachable from the source
 *   vertex in at most \c order steps are included. For \c IGRAPH_IN
 *   all vertices from which the source vertex is reachable in at most
 *   \c order steps are included. \c IGRAPH_ALL ignores the direction
 *   of the edges. This argument is ignored for undirected graphs.
 * \return Error code.
 * 
 * \sa \ref igraph_neighborhood_size() to calculate the size of the
 * neighborhood, \ref igraph_neighborhood_graphs() for creating
 * graphs from the neighborhoods.
 * 
 * Time complexity: O(n*d*o), n is the number of vertices for which
 * the calculation is performed, d is the average degree, o is the
 * order.
 */

int igraph_neighborhood(const igraph_t *graph,  // CSR MATRIX
												igraph_vector_ptr_t *res, // result vertices for some level of order
												igraph_vs_t vids, // input set
												igraph_integer_t order, //
												igraph_neimode_t mode
												)
{
  
  long int no_of_nodes=igraph_vcount(graph); // calculate number of nodes in the graph
  igraph_dqueue_t q; // double-ended queue (C++: std::deque)
	// igraph_dqueue_pop :=  Remove the head. (C++: pop_front)
	// igraph_dqueue_pop_back :=  Remove the tail. (C++: pop_back)
	// igraph_dqueue_push := Append an element to the end of the queue (head side) (C++: push_back) (there is no push_back)

	// igraph_dqueue_back := The last element in the queue (C++: back)
	// igraph_dqueue_head := The first element in the queue (C++: front)
	// igraph_dqueue_size := number of elements in the queue (C++: size)
	
  igraph_vit_t vit; // vertex iterator
  long int i, j;
  long int *added; // added vertices
  igraph_vector_t neis; //neighbors of a vertex
  igraph_vector_t tmp; // ?temporary help vector
  igraph_vector_t *newv; // ?new vertices

	// check order > 0
  
  added=igraph_Calloc(no_of_nodes, long int); // allocate memory for added nodes with size of "number of all nodes in the graph" of type "int"
	
  if (added==0) {
    IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added); // ~ probably free memory for added later on
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100); // initialize dqueue
  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit)); // create vit-iterator for the graph and iterate over input set
  IGRAPH_FINALLY(igraph_vit_destroy, &vit); // destroy iterator eventually
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0); // initialize and free neighbors-vector
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0); // initialize and free temporary help vector
  IGRAPH_CHECK(igraph_vector_ptr_resize(res, IGRAPH_VIT_SIZE(vit))); // ??? get size of the iterator and resize the initial result vertices-vector.
  
  for (i=0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) { //iterate over input set. (for each vid do .... )
    long int node=IGRAPH_VIT_GET(vit); // get the current vertex, aka node;
    added[node]=i+1; // add node to added vertices-vector in location "node (vertex_id)" and set it to vertexnumber + 1
    igraph_vector_clear(&tmp); // clear the temporary help-vector
    IGRAPH_CHECK(igraph_vector_push_back(&tmp, node)); // add current vertex to the back of the temporary help-vector.
    if (order > 0) {
      igraph_dqueue_push(&q, node); // add vertex to the end of the dqueue if order > 0
      igraph_dqueue_push(&q, 0); // add 0 afterwards to the end of the dqueue (this number is the distance of the node to itself)
    }

    while (!igraph_dqueue_empty(&q)) { // as long as dqueue is not empty do ...
      long int actnode=(long int) igraph_dqueue_pop(&q); // remove the head (the node)
      long int actdist=(long int) igraph_dqueue_pop(&q); // remove the head (the distance/order of the node)
      long int n; // declare n := number of neighbors
      igraph_neighbors(graph, &neis, (igraph_integer_t) actnode, mode); // get all neighbors (mode == out) from the current vertex and store in neis
      n=igraph_vector_size(&neis); // calculate size of neighbors
      
      if (actdist<order-1) { // ask if the node-distance is smaller than the order (in first iteration actdist == 0 and order-1 == 0 for order==1)
				/* we add them to the q */
				for (j=0; j<n; j++) { // iterate over all neighbors
					long int nei=(long int) VECTOR(neis)[j]; // get current neighbor
					if (added[nei] != i+1) { // only add if the added vertices-vector in location "nei (vertex_id)" is not the vertexnumber + 1
																	 // in other words this should mean if the neighbor is not the current vertex in the input set or has not been added yet!
																	 // *added is kind of like the marker for "has been visited"
						added[nei]=i+1; // mark as visited
						IGRAPH_CHECK(igraph_dqueue_push(&q, nei)); // add the vertex to the dqueue
						IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1)); // add its distance to the queue
						IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei)); // add the vertex to the result 
					}
				}
      } else {
				/* we just count them but don't add them to q */ // just add them to the result, but not process them in the dqueue
				for (j=0; j<n; j++) { // iterate over all neighbors
					long int nei=(long int) VECTOR(neis)[j];
					if (added[nei] != i+1) {
						added[nei]=i+1;
						IGRAPH_CHECK(igraph_vector_push_back(&tmp, nei));
					}
				}
      }
    } /* while q not empty */

    newv=igraph_Calloc(1, igraph_vector_t); // allocate memory for newv of size "1" of type "vector"
		// check calloc
		
    IGRAPH_FINALLY(igraph_free, newv);
    IGRAPH_CHECK(igraph_vector_copy(newv, &tmp));
    VECTOR(*res)[i]=newv;
    IGRAPH_FINALLY_CLEAN(1);
  } // end for every vertex in the input set

  igraph_vector_destroy(&tmp);
  igraph_vector_destroy(&neis);
  igraph_vit_destroy(&vit);
  igraph_dqueue_destroy(&q);
  igraph_Free(added);
  IGRAPH_FINALLY_CLEAN(5);

  return 0;
}