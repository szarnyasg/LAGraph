//------------------------------------------------------------------------------
// LAGraph_reorder_vertices: reorder vertices according to their degree.
//
// Vertex label reordering (also known as permutation) might lead to better
// performance for some algorithms, such as triangle count.
//
// ------------------------------------------------------------------------------

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

// LAGraph_reorder_vertices: Contributed by Tim Davis, Texas A&M and
// Gabor Szarnyas, Budapest University of Technology and Economics

#include "LAGraph_internal.h"
#include "GB_msort_2.h"

#define LAGRAPH_FREE_ALL        \
    GrB_free (&C) ;             \
    GrB_free (&L) ;             \
    GrB_free (&T) ;             \
    GrB_free (&U) ;             \
    LAGRAPH_FREE (W0) ;         \
    LAGRAPH_FREE (W1) ;         \
    LAGRAPH_FREE (P) ;          \
    LAGRAPH_FREE (D) ;

// TODO:
// currently only pattern matrices are supported

GrB_Info LAGraph_reorder_vertices // reorder vertices according to their degree
(
    GrB_Matrix *Chandle,         // output matrix
    GrB_Index *mapping,          // mapping from old vertices to new vertices
                                 // (unused if NULL)
    GrB_Matrix A,                // input matrix
    bool ascending               // sort in ascending order
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Info info ;

    GrB_Type type;
    GrB_Index n ;
    GrB_Vector degrees_vec = NULL;
    GrB_Index degrees_nvals ;
    GrB_Index* degree ;
    GrB_Matrix C = NULL, L = NULL, U = NULL, T = NULL ;
    int64_t *P = NULL, *D = NULL, *W0 = NULL, *W1 = NULL ;

    LAGr_Matrix_type (&type, A)
    LAGr_Matrix_nrows (&n, A)

    if (type != GrB_BOOL)
    {
        LAGraph_pattern(&C, A, GrB_BOOL);
    }
    else
    {
        C = A;
    }

    GrB_Vector X = NULL, Dv = NULL ;
    GrB_Index *Di = NULL ;
    LAGr_Vector_new (&X, GrB_BOOL, n) ;
    LAGr_Vector_new (&Dv, GrB_INT64, n) ;
    LAGr_assign (X, NULL, NULL, 0, GrB_ALL, n, NULL) ;
    LAGr_assign (Dv, NULL, NULL, 0, GrB_ALL, n, NULL) ;
    LAGr_vxm (Dv, NULL, GrB_PLUS_INT64, GxB_PLUS_PAIR_INT64, X, A, NULL) ;
    GrB_free (&X) ;
    // GxB_print (A, 2) ;
    // GxB_print (D, 2) ;
    GrB_free (&X) ;
    GrB_Index n2, nvals2 ;
    LAGr_Vector_export (&Dv, &type, &n2, &nvals2, &Di, (void **) &degree, NULL) ;
    if (n != n2 || n != nvals2) { printf ("??\n") ; abort ( ) ; }
    LAGRAPH_FREE (Di) ;


    #define CHUNK (64*1024)
    int nthreads = LAGraph_get_nthreads ( ) ;
    nthreads = LAGRAPH_MIN (nthreads, n/CHUNK) ;
    nthreads = LAGRAPH_MAX (nthreads, 1) ;

    // allocate workspace
    P  = LAGraph_malloc (n, sizeof (int64_t)) ;
    D  = LAGraph_malloc (n, sizeof (int64_t)) ;
    W0 = LAGraph_malloc (n, sizeof (int64_t)) ;
    W1 = LAGraph_malloc (n, sizeof (int64_t)) ;
    if (P == NULL || D == NULL || W0 == NULL || W1 == NULL)
    {
        // out of memory
        LAGRAPH_FREE_ALL ;
        return (GrB_OUT_OF_MEMORY) ;
    }

    // construct the pair [D,P] to sort
    if (ascending)
    {
        // sort [D,P] in ascending order of degree, tie-breaking on P
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int64_t k = 0 ; k < n ; k++)
        {
            D [k] = degree [k] ;
            P [k] = k ;
        }
    }
    else
    {
        // sort [D,P] in descending order of degree, tie-breaking on P
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int64_t k = 0 ; k < n ; k++)
        {
            D [k] = -degree [k] ;
            P [k] = k ;
        }
    }

    GB_msort_2 (D, P, W0, W1, n, nthreads) ;

    // T = A (P,P)
    LAGr_Matrix_new (&T, type, n, n) ;
    LAGr_extract (T, NULL, NULL, C, P, n, P, n, NULL) ;

    mapping = P ;
    (*Chandle) = T ;

    return (GrB_SUCCESS) ;
}
