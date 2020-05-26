//------------------------------------------------------------------------------
// LAGraph_msbfs: multi-source breadth-first search
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2020 LAGraph Contributors.

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

// LAGraph_msbfs: multi-source breadth-first search algorithm

// Usage:


// References:

// Carl Yang, Aydin Buluc, and John D. Owens. 2018. Implementing Push-Pull
// Efficiently in GraphBLAS. In Proceedings of the 47th International
// Conference on Parallel Processing (ICPP 2018). ACM, New York, NY, USA,
// Article 89, 11 pages. DOI: https://doi.org/10.1145/3225058.3225122

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&q) ;         \
    GrB_free (&L) ;         \
    GrB_free (&t) ;         \
    GrB_free (&pi) ;        \
}

GrB_Info LAGraph_msbfs      // multi-source BFS
(
    GrB_Matrix *L_handle,   // L(t,i) is the BFS levels of node i in traversal t
    GrB_Matrix A,           // input graph, treated as if boolean in semiring
    GrB_Matrix AT,          // transpose of A (optional; push-only if NULL)
    GrB_Index* s,           // starting nodes of the BFS
    GrB_Index ns,           // number of starting nodes
    int64_t max_level,      // optional limit of # levels to search
    bool vsparse            // if true, v is expected to be very sparse
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Vector q = NULL ;           // nodes visited at each level
    GrB_Matrix L = NULL ;           // result matrix
    GrB_Vector t = NULL ;           // temporary vector
    GrB_Vector pi = NULL ;          // parent vector

    if (L_handle == NULL || (A == NULL && AT == NULL))
    {
        // required output argument is missing
        LAGRAPH_ERROR ("required arguments are NULL", GrB_NULL_POINTER) ;
    }

    (*L_handle) = NULL ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS ) \
        && ( GxB_IMPLEMENTATION >= GxB_VERSION (3,2,0) )
    GrB_Descriptor desc_s  = GrB_DESC_S ;
    GrB_Descriptor desc_sc = GrB_DESC_SC ;
    GrB_Descriptor desc_rc = GrB_DESC_RC ;
    GrB_Descriptor desc_r  = GrB_DESC_R ;
    #else
    GrB_Descriptor desc_s  = NULL ;
    GrB_Descriptor desc_sc = LAGraph_desc_ooco ;
    GrB_Descriptor desc_rc = LAGraph_desc_oocr ;
    GrB_Descriptor desc_r  = LAGraph_desc_ooor ;
    #endif

    bool use_vxm_with_A ;
    GrB_Index nrows, ncols, nvalA, ignore, nvals ;
    if (A == NULL)
    {
        // only AT is provided
        LAGr_Matrix_ncols (&nrows, AT) ;
        LAGr_Matrix_nrows (&ncols, AT) ;
        LAGr_Matrix_nvals (&nvalA, AT) ;
        use_vxm_with_A = false ;
    }
    else
    {
        // A is provided.  AT may or may not be provided
        LAGr_Matrix_nrows (&nrows, A) ;
        LAGr_Matrix_ncols (&ncols, A) ;
        LAGr_Matrix_nvals (&nvalA, A) ;
        use_vxm_with_A = true ;
    }

    // push/pull requires both A and AT
    bool push_pull = (A != NULL && AT != NULL) ;

    if (nrows != ncols)
    {
        // A must be square
        LAGRAPH_ERROR ("A must be square", GrB_NULL_POINTER) ;
    }

    //--------------------------------------------------------------------------
    // check the format of A and AT
    //--------------------------------------------------------------------------

    // TODO



    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*L_handle) = L ;       // return result
    L = NULL ;              // set to NULL so LAGRAPH_FREE_ALL doesn't free it
    LAGRAPH_FREE_ALL ;      // free all workspace (except for result v)
    return (GrB_SUCCESS) ;
}

