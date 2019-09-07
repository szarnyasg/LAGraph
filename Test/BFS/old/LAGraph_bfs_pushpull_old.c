//------------------------------------------------------------------------------
// LAGraph_bfs_pushpull_old:  push-pull BFS
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

// Contributed by Tim Davis, Texas A&M

// hacked, older version, of LAGraph_bfs_pushpull.  Delete this.

#include "bfs_test.h"

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&v) ;         \
    GrB_free (&q) ;         \
    GrB_free (&desc) ;      \
}

#ifdef MY_BOOL
// just to test compile-time user-defined objects in SuiteSparse/GraphBLAS:
#define Boolean My_LOR_LAND
#else
#define Boolean LAGraph_LOR_LAND_BOOL
#endif

GrB_Info LAGraph_bfs_pushpull_old
(
    GrB_Vector *v_output,   // v [i] is the BFS level of node i in the graph
    const GrB_Matrix A,     // input graph, treated as if boolean in semiring
    const GrB_Matrix AT,    // transpose of A
    GrB_Index s,            // starting node of the BFS
    int32_t max_level       // max # of levels to search (<0: nothing,
                            // 1: just the source, 2: source and neighbors, etc)
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    #ifdef MY_BOOL
    fprintf (stderr, "\npush-pull (user-defined):\n") ;
    #endif

    GrB_Info info ;

    GrB_Vector q = NULL ;           // nodes visited at each level
    GrB_Vector v = NULL ;           // result vector
    GrB_Descriptor desc = NULL ;    // descriptor for vxm

    if (v_output == NULL)
    {
        // required output argument is missing
        return (GrB_NULL_POINTER) ;
    }

    (*v_output) = NULL ;

    GrB_Index nrows, ncols, nvalsA ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A)) ;
    LAGRAPH_OK (GrB_Matrix_nvals (&nvalsA, A)) ;
    if (nrows != ncols)
    {
        // A must be square
        return (GrB_INVALID_VALUE) ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index n = nrows ;
    max_level = LAGRAPH_MIN (n, max_level) ;

    // create an empty vector v.  Assume int32 is sufficient
    LAGRAPH_OK (GrB_Vector_new (&v, GrB_INT32, n)) ;

    // make it dense
    // if (max_level > 2)
    {
        for (GrB_Index i = 0 ; i < n ; i++)
        {
            LAGRAPH_OK (GrB_Vector_setElement (v, (int32_t) 0, i)) ;
        }
    }

    // create a boolean vector q, and set q[s] to true
    LAGRAPH_OK (GrB_Vector_new (&q, GrB_BOOL, n)) ;
    LAGRAPH_OK (GrB_Vector_setElement (q, true, s)) ;

            // hack
//          GrB_Descriptor_set (LAGraph_desc_oocr, GxB_AxB_METHOD ,
//              GxB_AxB_GUSTAVSON) ;

    //--------------------------------------------------------------------------
    // BFS traversal and label the nodes
    //--------------------------------------------------------------------------

    char filename [1000] ;
    sprintf (filename, "results/bfs_pushpull_n%d.m", (int) n) ;
    FILE *f = fopen (filename, "w") ;
    double tic [2], t ;
    fprintf (f, "\nfunction [a m s q] = bfs_pushpull_n%d\n", (int) n) ;

    int64_t nvisited = 0 ;      // # nodes visited so far
    int64_t nvstart = 0 ;

    // average node degree
    double d = (n == 0) ? 0 : (((double) nvalsA) / ((double) n)) ;

    bool successor = true ; // true when some successor found
    for (int32_t level = 1 ; successor && level <= max_level ; level++)
    {

        // TODO  let v start sparse.  when it reaches
        // nvals(v) > alpha*sqrt(n), convert to dense.
        // find a good default alpha.  Provide alpha on input:
        //  alpha=0:  use a dense v to start with
        //  alhpa<0:  automatic selection
        //  alpha>0:  use this particular value for alpha.

        GrB_Index nq ;
        LAGRAPH_OK (GrB_Vector_nvals (&nq, q)) ;
        fprintf (f, "q(%5g) = %10g ;", (double) level, (double) nq) ;

        // v<q> = level, using vector assign with q as the mask
        LAGraph_tic (tic) ;
        LAGRAPH_OK (GrB_assign (v, q, NULL, level, GrB_ALL, n, NULL)) ;
        t = LAGraph_toc (tic) ;
        fprintf (f, "a(%5g) = %12.6e ; ", (double) level, t) ;

        nvisited += nq ;
        // fprintf (stderr, "level %g nvisited %g nq %g\n",
        //     (double) level, (double) nvisited, (double) nq) ;

        if (nvisited == n)
        {
            // fprintf (stderr, "break HERE\n\n") ;
            fprintf (f, "m(%5g) = 0 ; ", (double) level) ;
            fprintf (f, "s(%5g) = 0 ;\n", (double) level) ;
            break ;
        }

        // GxB_fprint (AT, 3, stderr) ;
        // fprintf (stderr, "v is now:\n") ;
        // GxB_fprint (v, 3, stderr) ;
        // GxB_fprint (LAGraph_desc_oocr, 3, stderr) ;
        // GxB_fprint (LAGraph_LOR_LAND_BOOL, 3, stderr) ;

        // LAGRAPH_OK (GrB_Vector_nvals (&nvals, q)) ;

            double pushwork = d * ((double) nq) ;
            double expected = (double) n / (double) (nvisited-nvstart+1) ;
            double per_dot = LAGRAPH_MIN (d, expected) ;
            double binarysearch = (3 * (1 + log2 ((double) nq))) ;
            double pullwork = (n-nvisited) * per_dot * binarysearch ;

            // fprintf (stderr, "\nlevel %g\n", (double) level) ;
            // fprintf (stderr, "pushwork %g\n", pushwork) ;
            // fprintf (stderr, "expected %g\n", expected) ;
            // fprintf (stderr, "per_dot  %g\n", per_dot) ;
            // fprintf (stderr, "binary   %g\n", binarysearch) ;
            // fprintf (stderr, "pullwork %g\n", pullwork) ;

        // q<!v> = q ||.&& A ; finds all the unvisited
        // successors from current q, using !v as the mask
        LAGraph_tic (tic) ;
        if (pushwork < pullwork)
        {

           // fprintf (stderr, "%g %g : did push %g\n",
           //     pushwork, pullwork, (double) level) ;

            // push, using saxpy operations
            // TODO in Version 2.2.3 of SuiteSparse:GraphBLAS, a dense mask
            // vector will be slow (the entire vector gets scattered, but
            // not all is used...).  This is fixed in 2.3.0 (not released)

           // fprintf (stderr, "%g %g : did push %g\n",
           // pushwork, pullwork, (double) level) ;

            LAGRAPH_OK (GrB_vxm (q, v, NULL, LAGraph_LOR_LAND_BOOL, q, A,
                LAGraph_desc_oocr)) ;
        }
        else
        {
            // pull, using dot products
            // TODO this needs early-exit (not on v2.2.3 of SuiteSparse
            // but will be in the next release)
            // fprintf (stderr, "level %g: --- pull (dot):\n",
            //    (double) level) ;

           // fprintf (stderr, "%g %g : did pull %g\n",
           // pushwork, pullwork, (double) level) ;

            // q<!v> = AT*q
            LAGRAPH_OK (GrB_mxv (q, v, NULL, LAGraph_LOR_LAND_BOOL, AT, q,
                LAGraph_desc_oocr)) ;

            #if 0
            // q<!v> = AT*v ; gives the same result.  If v(j) != 0, then
            // node i has been seen.  AT(i,j)*v(j) is true if the edge (j,i)
            // is present and v(j) is nonzero.  This sets y(i) to true
            // where y = AT*v.  However, the value y(i) is assigned to q(i)
            // only if v(i) was initially zero.
            LAGRAPH_OK (GrB_mxv (q, v, NULL, LAGraph_LOR_LAND_BOOL, AT, v,
                LAGraph_desc_oocr)) ;

            // make q sparse
            LAGRAPH_OK (GrB_assign (q, q, NULL, q, GrB_ALL, n,
                LAGraph_desc_ooor)) ;
            #endif
        }

        t = LAGraph_toc (tic) ;
        fprintf (f, "m (%5g) = %12.6e ; ", (double) level, t) ;

        // GxB_fprint (q, 3, stderr) ;

        // dump q to see how it was computed:
        // GxB_fprint (q, GxB_COMPLETE, stderr) ;

        // note that if A has no explicit zeros, then this works faster:
        // GrB_Vector_nvals (&nvals, q) ; successor = (nvals > 0) ;

        // successor = ||(q)
        LAGraph_tic (tic) ;
        LAGRAPH_OK (GrB_reduce (&successor, NULL, LAGraph_LOR_MONOID, q, NULL));
        t = LAGraph_toc (tic) ;
        fprintf (f, "s (%5g) = %12.6e ;\n", (double) level, t) ;
    }

    fclose (f) ;

    //--------------------------------------------------------------------------
    // make v sparse
    //--------------------------------------------------------------------------

    GrB_Index nvals ;
    LAGRAPH_OK (GrB_Vector_nvals (&nvals, v)) ;

    if (nvals < n)
    {
        // v<v> = v ; clearing v before assigning it back
        LAGRAPH_OK (GrB_assign (v, v, NULL, v, GrB_ALL, n, LAGraph_desc_ooor)) ;
    }

    // finish the work
    LAGRAPH_OK (GrB_Vector_nvals (&nvals, v)) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*v_output) = v ;       // return result
    v = NULL ;              // set to NULL so LAGRAPH_FREE_ALL doesn't free it
    LAGRAPH_FREE_ALL ;      // free all workspace (except for result v)
    return (GrB_SUCCESS) ;
}

