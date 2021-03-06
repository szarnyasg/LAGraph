//------------------------------------------------------------------------------
// LAGraph_BF_full2.c: Bellman-Ford single-source shortest paths, returns tree,
// while diagonal of input matrix A needs not to be explicit 0, using the
// frontier idea from Roi Lipman
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// LAGraph_BF_full2: Bellman-Ford single source shortest paths, returning both
// the path lengths and the shortest-path tree.  contributed by Jinhao Chen and
// Tim Davis, Texas A&M.

// LAGraph_BF_full2 performs a Bellman-Ford to find out shortest path, parent
// nodes along the path and the hops (number of edges) in the path from given
// source vertex s in the range of [0, n) on graph given as matrix A with size
// n*n. The sparse matrix A has entry A(i, j) if there is an edge from vertex i
// to vertex j with weight w, then A(i, j) = w.

// TODO: think about the return values

// LAGraph_BF_full2 returns GrB_SUCCESS if it succeeds.  In this case, there
// are no negative-weight cycles in the graph, and d, pi, and h are returned.
// The vector d has d(k) as the shortest distance from s to k. pi(k) = p+1,
// where p is the parent node of k-th node in the shortest path. In particular,
// pi(s) = 0. h(k) = hop(s, k), the number of edges from s to k in the shortest
// path.

// If the graph has a negative-weight cycle, GrB_NO_VALUE is returned, and the
// GrB_Vectors d(k), pi(k) and h(k)  (i.e., *pd_output, *ppi_output and
// *ph_output respectively) will be NULL when negative-weight cycle detected.

// Otherwise, other errors such as GrB_OUT_OF_MEMORY, GrB_INVALID_OBJECT, and
// so on, can be returned, if these errors are found by the underlying
// GrB_* functions.

//------------------------------------------------------------------------------
#include "BF_test.h"

#define LAGRAPH_FREE_WORK              \
{                                      \
    GrB_free(&d);                      \
    GrB_free(&dtmp);                   \
    GrB_free(&dfrontier);              \
    GrB_free(&Atmp);                   \
    GrB_free(&BF_Tuple3);              \
    GrB_free(&BF_lMIN_Tuple3);         \
    GrB_free(&BF_PLUSrhs_Tuple3);      \
    GrB_free(&BF_EQ_Tuple3);           \
    GrB_free(&BF_lMIN_Tuple3_Monoid);  \
    GrB_free(&BF_lMIN_PLUSrhs_Tuple3); \
    LAGRAPH_FREE (I);                  \
    LAGRAPH_FREE (J);                  \
    LAGRAPH_FREE (w);                  \
    LAGRAPH_FREE (W);                  \
    LAGRAPH_FREE (h);                  \
    LAGRAPH_FREE (pi);                 \
} 
 
#define LAGRAPH_FREE_ALL               \
{                                      \
    LAGRAPH_FREE_WORK                  \
    GrB_free (pd_output);              \
    GrB_free (ppi_output);             \
    GrB_free (ph_output);              \
}


//------------------------------------------------------------------------------
// data type for each entry of the adjacent matrix A and "distance" vector d;
// <INFINITY,INFINITY,INFINITY> corresponds to nonexistence of a path, and
// the value  <0, 0, NULL> corresponds to a path from a vertex to itself
//------------------------------------------------------------------------------
typedef struct
{
    double w;    // w  corresponds to a path weight.
    GrB_Index h; // h  corresponds to a path size or number of hops.
    GrB_Index pi;// pi corresponds to the penultimate vertex along a path.
                 // vertex indexed as 1, 2, 3, ... , V, and pi = 0 (as nil)
                 // for u=v, and pi = UINT64_MAX (as inf) for (u,v) not in E
}
BF_Tuple3_struct;

//------------------------------------------------------------------------------
// 2 binary functions, z=f(x,y), where Tuple3xTuple3 -> Tuple3
//------------------------------------------------------------------------------
void BF_lMIN2
(
    BF_Tuple3_struct *z,
    const BF_Tuple3_struct *x,
    const BF_Tuple3_struct *y
)
{
    if (x->w < y->w
        || (x->w == y->w && x->h < y->h)
        || (x->w == y->w && x->h == y->h && x->pi < y->pi))
    {
        if (z != x) { *z = *x; }
    }
    else
    {
        *z = *y;
    }
}

void BF_PLUSrhs2
(
    BF_Tuple3_struct *z,
    const BF_Tuple3_struct *x,
    const BF_Tuple3_struct *y
)
{
    z->w = x->w + y->w;
    z->h = x->h + y->h;
    if (x->pi != UINT64_MAX && y->pi != 0)
    {
        z->pi = y->pi;
    }
    else
    {
        z->pi = x->pi;
    }
}

void BF_EQ
(
    bool *z,
    const BF_Tuple3_struct *x,
    const BF_Tuple3_struct *y
)
{
    if (x->w == y->w && x->h == y->h && x->pi == y->pi)
    {
        *z = true;
    }
    else
    {
        *z = false;
    }
}

// Given a n-by-n adjacency matrix A and a source vertex s.
// If there is no negative-weight cycle reachable from s, return the distances
// of shortest paths from s and parents along the paths as vector d. Otherwise,
// returns d=NULL if there is a negtive-weight cycle.
// pd_output is pointer to a GrB_Vector, where the i-th entry is d(s,i), the
//   sum of edges length in the shortest path
// ppi_output is pointer to a GrB_Vector, where the i-th entry is pi(i), the
//   parent of i-th vertex in the shortest path
// ph_output is pointer to a GrB_Vector, where the i-th entry is h(s,i), the
//   number of edges from s to i in the shortest path
// A has weights on corresponding entries of edges
// s is given index for source vertex
GrB_Info LAGraph_BF_full2
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    GrB_Vector *ppi_output,     //the pointer to the vector of parent
    GrB_Vector *ph_output,       //the pointer to the vector of hops
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Index s           //given index of the source
)
{
    GrB_Info info;
    // tmp vector to store distance vector after n (i.e., V) loops
    GrB_Vector d = NULL, dtmp = NULL, dfrontier = NULL;
    GrB_Matrix Atmp = NULL;
    GrB_Type BF_Tuple3;

    GrB_BinaryOp BF_lMIN_Tuple3;
    GrB_BinaryOp BF_PLUSrhs_Tuple3;
    GrB_BinaryOp BF_EQ_Tuple3;

    GrB_Monoid BF_lMIN_Tuple3_Monoid;
    GrB_Semiring BF_lMIN_PLUSrhs_Tuple3;

    GrB_Index nrows, ncols, n, nz;  // n = # of row/col, nz = # of nnz in graph
    GrB_Index *I = NULL, *J = NULL; // for col/row indices of entries from A
    GrB_Index *h = NULL, *pi = NULL;
    double *w = NULL;
    BF_Tuple3_struct *W = NULL;

    if (A == NULL || pd_output == NULL ||
        ppi_output == NULL || ph_output == NULL)
    {
        // required argument is missing
        LAGRAPH_ERROR ("required arguments are NULL", GrB_NULL_POINTER) ;
    }

    *pd_output  = NULL;
    *ppi_output = NULL;
    *ph_output  = NULL;
    LAGr_Matrix_nrows (&nrows, A) ;
    LAGr_Matrix_ncols (&ncols, A) ;
    LAGr_Matrix_nvals (&nz, A);
    if (nrows != ncols)
    {
        // A must be square
        LAGRAPH_ERROR ("A must be square", GrB_INVALID_VALUE) ;
    }
    n = nrows;

    if (s >= n || s < 0)
    {
        LAGRAPH_ERROR ("invalid value for source vertex s", GrB_INVALID_VALUE);
    }
    //--------------------------------------------------------------------------
    // create all GrB_Type GrB_BinaryOp GrB_Monoid and GrB_Semiring
    //--------------------------------------------------------------------------
    // GrB_Type
    LAGr_Type_new(&BF_Tuple3, sizeof(BF_Tuple3_struct));

    // GrB_BinaryOp
    LAGr_BinaryOp_new(&BF_EQ_Tuple3,
        (LAGraph_binary_function) (&BF_EQ), GrB_BOOL, BF_Tuple3, BF_Tuple3);
    LAGr_BinaryOp_new(&BF_lMIN_Tuple3,
        (LAGraph_binary_function) (&BF_lMIN2),
        BF_Tuple3, BF_Tuple3, BF_Tuple3);
    LAGr_BinaryOp_new(&BF_PLUSrhs_Tuple3,
        (LAGraph_binary_function)(&BF_PLUSrhs2),
        BF_Tuple3, BF_Tuple3, BF_Tuple3);

    // GrB_Monoid
    BF_Tuple3_struct BF_identity = (BF_Tuple3_struct) { .w = INFINITY,
        .h = UINT64_MAX, .pi = UINT64_MAX };
    LAGRAPH_OK(GrB_Monoid_new_UDT(&BF_lMIN_Tuple3_Monoid, BF_lMIN_Tuple3,
        &BF_identity));

    //GrB_Semiring
    LAGr_Semiring_new(&BF_lMIN_PLUSrhs_Tuple3,
        BF_lMIN_Tuple3_Monoid, BF_PLUSrhs_Tuple3);

    //--------------------------------------------------------------------------
    // allocate arrays used for tuplets
    //--------------------------------------------------------------------------
    I = LAGraph_malloc (nz, sizeof(GrB_Index)) ;
    J = LAGraph_malloc (nz, sizeof(GrB_Index)) ;
    w = LAGraph_malloc (nz, sizeof(double)) ;
    W = LAGraph_malloc (nz, sizeof(BF_Tuple3_struct)) ;
    if (I == NULL || J == NULL || w == NULL || W == NULL)
    {
        LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
    }

    //--------------------------------------------------------------------------
    // create matrix Atmp based on A, while its entries become BF_Tuple3 type
    //--------------------------------------------------------------------------
    LAGRAPH_OK(GrB_Matrix_extractTuples_FP64(I, J, w, &nz, A));
    int nthreads = LAGraph_get_nthreads ( ) ;
    printf ("nthreads %d\n", nthreads) ;
    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (GrB_Index k = 0; k < nz; k++)
    {
        if (w[k] == 0)             //diagonal entries
        {
            W[k] = (BF_Tuple3_struct) { .w = 0, .h = 0, .pi = 0 };
        }
        else
        {
            W[k] = (BF_Tuple3_struct) { .w = w[k], .h = 1, .pi = I[k] + 1 };
        }
    }
    LAGr_Matrix_new(&Atmp, BF_Tuple3, n, n);
    LAGRAPH_OK(GrB_Matrix_build_UDT(Atmp, I, J, W, nz, BF_lMIN_Tuple3));
    LAGRAPH_FREE (I);
    LAGRAPH_FREE (J);
    LAGRAPH_FREE (W);
    LAGRAPH_FREE (w);

    //--------------------------------------------------------------------------
    // create and initialize "distance" vector d
    //--------------------------------------------------------------------------
    LAGr_Vector_new(&d, BF_Tuple3, n);
    // initial distance from s to itself
    BF_Tuple3_struct d0 = (BF_Tuple3_struct) { .w = 0, .h = 0, .pi = 0 };
    LAGRAPH_OK(GrB_Vector_setElement_UDT(d, &d0, s));

    //--------------------------------------------------------------------------
    // start the Bellman Ford process
    //--------------------------------------------------------------------------
    // copy d to dtmp in order to create a same size of vector
    LAGr_Vector_dup(&dtmp, d);
    LAGr_Vector_dup(&dfrontier, d);
    bool same= false;          // variable indicating if d == dtmp
    int64_t iter = 0;          // number of iterations

    // terminate when no new path is found or more than V-1 loops
    while (!same && iter < n - 1)
    {
        // execute semiring on d and A, and save the result to dtmp
        LAGr_vxm(dfrontier, GrB_NULL, GrB_NULL,
            BF_lMIN_PLUSrhs_Tuple3, dfrontier, Atmp, GrB_NULL);

        // dtmp[i] = min(d[i], dfrontier[i]).
        GrB_eWiseAdd_Vector_BinaryOp(dtmp, GrB_NULL, GrB_NULL, BF_lMIN_Tuple3,
            d, dfrontier, GrB_NULL);
        
        LAGRAPH_OK (LAGraph_Vector_isequal(&same, dtmp, d, BF_EQ_Tuple3));
        if (!same)
        {
            GrB_Vector ttmp = dtmp;
            dtmp = d;
            d = ttmp;
        }
        iter ++;
    }

    // check for negative-weight cycle only when there was a new path in the
    // last loop, otherwise, there can't be a negative-weight cycle.
    if (!same)
    {
        // execute semiring again to check for negative-weight cycle
        LAGr_vxm(dfrontier, GrB_NULL, GrB_NULL,
            BF_lMIN_PLUSrhs_Tuple3, dfrontier, Atmp, GrB_NULL);

        // dtmp[i] = min(d[i], dfrontier[i]).
        GrB_eWiseAdd_Vector_BinaryOp(dtmp, GrB_NULL, GrB_NULL, BF_lMIN_Tuple3,
            d, dfrontier, GrB_NULL);

        // if d != dtmp, then there is a negative-weight cycle in the graph
        LAGRAPH_OK (LAGraph_Vector_isequal(&same, dtmp, d, BF_EQ_Tuple3));
        if (!same)
        {
            // printf("A negative-weight cycle found. \n");
            LAGRAPH_FREE_ALL;
            return (GrB_NO_VALUE) ;
        }
    }

    //--------------------------------------------------------------------------
    // extract tuple from "distance" vector d and create GrB_Vectors for output
    //--------------------------------------------------------------------------
    I = LAGraph_malloc (n, sizeof(GrB_Index)) ;
    W = LAGraph_malloc (n, sizeof(BF_Tuple3_struct)) ;
    w = LAGraph_malloc (n, sizeof(double)) ;
    h  = LAGraph_malloc (n, sizeof(GrB_Index)) ;
    pi = LAGraph_malloc (n, sizeof(GrB_Index)) ;
    if (I == NULL || W == NULL || w == NULL || h == NULL || pi == NULL)
    {
        LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
    }

    LAGRAPH_OK(GrB_Vector_extractTuples_UDT (I, (void *) W, &nz, d));

    for (GrB_Index k = 0; k < nz; k++)
    {
        w [k] = W[k].w ;
        h [k] = W[k].h ;
        pi[k] = W[k].pi;
    }
    LAGr_Vector_new(pd_output,  GrB_FP64,   n);
    LAGr_Vector_new(ppi_output, GrB_UINT64, n);
    LAGr_Vector_new(ph_output,  GrB_UINT64, n);
    LAGr_Vector_build (*pd_output , I, w , nz,GrB_MIN_FP64  );
    LAGr_Vector_build (*ppi_output, I, pi, nz,GrB_MIN_UINT64);
    LAGr_Vector_build (*ph_output , I, h , nz,GrB_MIN_UINT64);
    LAGRAPH_FREE_WORK;
    return (GrB_SUCCESS) ;
}
