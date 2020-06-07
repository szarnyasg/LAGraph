//------------------------------------------------------------------------------
// LAGraph/Test/MSBFS/msbfstest.c: test program for LAGraph_msbfs
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

// Contributed by Tim Davis, Texas A&M and Gabor Szarnyas, BME

// Usage: msbfstest ...
//

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL                            \
{                                                   \
    GrB_free (&A) ;                                 \
}

uint64_t extract_v(const GrB_Vector v, const GrB_Index i) {
    uint64_t x;
    GrB_Info info = GrB_Vector_extractElement(&x, v, i);
    if (info == GrB_NO_VALUE) {
        return 0xAAAAAAAAAAAAAAAA;
    }
    return x;
}

void print_bit_matrices(const GrB_Vector frontier, const GrB_Vector next, const GrB_Vector seen) {
    GrB_Index size;
    GrB_Vector_size(&size, frontier);
    printf("          frontier             next               seen\n");
    for (GrB_Index i = 0; i < size; i++) {
        printf("%4ld:", i);
        printf(
            " %016lx   %016lx   %016lx",
            extract_v(frontier, i), extract_v(next, i), extract_v(seen, i)
        );
        printf("\n");
    }
    printf("\n");
}

int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Matrix A = NULL;
    GrB_Vector frontier = NULL, next = NULL, seen = NULL, not_seen = NULL, filtered = NULL;

    LAGraph_init ( ) ;

    GrB_Semiring LAGr_BOR_FIRST = NULL, LAGr_BOR_SECOND = NULL ;
    LAGRAPH_OK (GrB_Semiring_new(&LAGr_BOR_FIRST, GxB_BOR_UINT64_MONOID, GrB_FIRST_UINT64))
    LAGRAPH_OK (GrB_Semiring_new(&LAGr_BOR_SECOND, GxB_BOR_UINT64_MONOID, GrB_SECOND_UINT64))

    // create the input matrix
    const GrB_Index n = 6;

    LAGr_Matrix_new(&A, GrB_BOOL, n, n)
    const bool val = true;

    // upper triangle
    LAGr_Matrix_setElement(A, val, 0, 2)
    LAGr_Matrix_setElement(A, val, 0, 3)
    LAGr_Matrix_setElement(A, val, 1, 2)
    LAGr_Matrix_setElement(A, val, 1, 3)
    LAGr_Matrix_setElement(A, val, 2, 4)
    LAGr_Matrix_setElement(A, val, 3, 5)
//    LAGr_Matrix_setElement(A, val, 6, 7)
    // lower triangle
    GrB_eWiseAdd(A, NULL, NULL, GrB_PLUS_UINT64, A, A, LAGraph_desc_otoo);

    LAGr_Vector_new(&frontier, GrB_UINT64, n)
    LAGr_Vector_new(&next, GrB_UINT64, n)
    LAGr_Vector_new(&seen, GrB_UINT64, n)
    LAGr_Vector_new(&not_seen, GrB_UINT64, n)
    LAGr_Vector_new(&filtered, GrB_UINT64, n)

    // wavefronts meeting
    GxB_Scalar meeting = NULL ;
    // Support scalar for GxB_select
    LAGRAPH_OK (GxB_Scalar_new (&meeting, GrB_UINT64))
    LAGRAPH_OK (GxB_Scalar_setElement (meeting, 0xC000000000000000L))


    // initialize frontier and seen matrices: to compute bidirectional search, start off with two searches
    GrB_Index* I = LAGraph_malloc(n, sizeof(GrB_Index));
    uint64_t*  X = LAGraph_malloc(n, sizeof(uint64_t));
    I[0] = 4; I[1] = 5; // even length path
//    I[0] = 4; I[1] = 3; // odd length path
    X[0] = 0x8000000000000000L;
    X[1] = 0x4000000000000000L;
    GrB_Vector_build(frontier, I, X, 2, GrB_BOR_UINT64);

    LAGr_Vector_dup(&seen, frontier)

    // initialize

    // traversal
    for (GrB_Index level = 1; level < n; level++) {
        printf("========================= Level %2ld =========================\n\n", level);

        // next = frontier * A
        bool push = true; // TODO: add heuristic
        if (push) {
            LAGr_vxm(next, NULL, NULL, LAGr_BOR_FIRST, frontier, A, NULL)
        } else {
            LAGr_mxv(next, NULL, NULL, LAGr_BOR_SECOND, A, frontier, NULL)
        }

        LAGr_apply(not_seen, NULL, NULL, GrB_BNOT_UINT64, seen, NULL)
//        print_bit_matrices(frontier, next, not_seen);
        // next = next & ~seen
        LAGr_eWiseAdd(next, next, NULL, GrB_BAND_UINT64, next, not_seen, NULL)

        GrB_Index next_nvals;
        LAGr_Vector_nvals(&next_nvals, next)
        if (next == 0) {
            printf("no new vertices found\n");
            break;
        }

        print_bit_matrices(frontier, next, seen);

        // even length path
        GrB_Index filtered_nvals;
        LAGr_select(filtered, NULL, NULL, GxB_EQ_THUNK, next, meeting, NULL)
        LAGr_Vector_nvals(&filtered_nvals, filtered)
        if (filtered_nvals > 0) { printf("found even length path: %ld\n", 2*level); break; }

        // seen = seen | next
        LAGr_eWiseAdd(seen, NULL, NULL, GrB_BOR_UINT64, seen, next, NULL)

        // odd length path
        LAGr_select(filtered, NULL, NULL, GxB_EQ_THUNK, seen, meeting, NULL)
        LAGr_Vector_nvals(&filtered_nvals, filtered)
        if (filtered_nvals > 0) { printf("found odd length path: %ld\n", 2*level-1); break; }

        // frontier = next
        LAGr_Vector_dup(&frontier, next)
    }

    LAGRAPH_FREE_ALL ;
    LAGraph_finalize ( ) ;
}
