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

uint64_t extract(const GrB_Matrix A, const GrB_Index i, const GrB_Index j) {
    uint64_t x;
    GrB_Info info = GrB_Matrix_extractElement(&x, A, i, j);
    if (info == GrB_NO_VALUE) {
        return 0;
    }
    return x;
}

uint64_t extract_v(const GrB_Vector v, const GrB_Index i) {
    uint64_t x;
    GrB_Info info = GrB_Vector_extractElement(&x, v, i);
    if (info == GrB_NO_VALUE) {
        return 0;
    }
    return x;
}

void print_bit_matrix(const GrB_Matrix A) {
    GrB_Index nrows, ncols;
    GrB_Matrix_nrows(&nrows, A);
    GrB_Matrix_ncols(&ncols, A);
    for (GrB_Index i = 0; i < nrows; i++) {
        printf("%4ld:", i);
        for (GrB_Index j = 0; j < ncols; j++) {
            printf(" %016lx", extract(A, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void print_bit_matrices(const GrB_Matrix frontier, const GrB_Matrix next, const GrB_Matrix seen, const GrB_Vector popcount, const GrB_Vector sp) {
    GrB_Index nrows, ncols;
    GrB_Matrix_nrows(&nrows, frontier);
    GrB_Matrix_ncols(&ncols, frontier);
    printf("          frontier             next               seen         popcount    sp\n");
    for (GrB_Index i = 0; i < nrows; i++) {
        printf("%4ld:", i);
        for (GrB_Index j = 0; j < ncols; j++) {
            printf(" %016lx   %016lx   %016lx   %8ld   %3ld", extract(frontier, i, j), extract(next, i, j), extract(seen, i, j), extract_v(popcount, i), extract_v(sp, i));
        }
        printf("\n");
    }
    printf("\n");
}

GrB_Info create_diagonal_bit_matrix(GrB_Matrix D) {
    GrB_Info info;

    GrB_Index n;
    GrB_Matrix_nrows(&n, D);

//    I = 0, 1, ..., n
//    J = 0, 0, ..., 0 [64], 1, 1, ..., 1 [64], ..., ceil(n/64)
//    X = repeat {b100..., b010..., b001..., ⋯, b...001} until we have n elements
    GrB_Index* I = LAGraph_malloc(n, sizeof(GrB_Index));
    GrB_Index* J = LAGraph_malloc(n, sizeof(GrB_Index));
    uint64_t*  X = LAGraph_malloc(n, sizeof(uint64_t));

    // TODO: parallelize this
    for (GrB_Index k = 0; k < n; k++) {
        I[k] = k;
        J[k] = k/64;
        X[k] = 0x8000000000000000L >> (k%64);
    }
    GrB_Matrix_build(D, I, J, X, n, GrB_BOR_UINT64);

    return info;
}

GrB_Info create_ones_vector(GrB_Vector v) {
    GrB_Info info;

    GrB_Index n;
    GrB_Vector_size(&n, v);

//    I = 0, 1, ..., n
//    X = 1, 1, ..., 1
    GrB_Index* I = LAGraph_malloc(n, sizeof(GrB_Index));
    uint64_t*  X = LAGraph_malloc(n, sizeof(uint64_t));

    for (GrB_Index k = 0; k < n; k++) {
        I[k] = k;
        X[k] = 1L;
    }
    GrB_Vector_build(v, I, X, n, GrB_PLUS_UINT64);

    return info;
}


uint64_t popcount(uint64_t x)
{
    int c = 0;
    for (; x != 0; x >>= 1)
        if (x & 1)
            c++;
    return c;
}

void fun_sum_popcount (void *z, const void *x)
{
    (*((uint64_t *) z))  = popcount(* ((uint64_t *) x));
}


int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Matrix A = NULL, frontier = NULL, next = NULL, seen = NULL, not_seen = NULL ;

    GrB_Matrix PopCount;
    GrB_Vector popcount;
    uint64_t total_popcount;

    GrB_Vector ones, level_v, sp;

    LAGraph_init ( ) ;

    // initializing unary operator for popcount
    GrB_UnaryOp op_popcount = NULL ;
    LAGRAPH_OK (GrB_UnaryOp_new(&op_popcount, fun_sum_popcount, GrB_UINT64, GrB_UINT64))

    // create the input matrix
    const GrB_Index n = 6;
    const GrB_Index bit_matrix_ncols = (n+63)/64;

    LAGr_Matrix_new(&A, GrB_UINT64, n, n)
    const uint64_t val = (uint64_t) -1; // 0xffff ffff ffff ffff

    // upper triangle
    LAGr_Matrix_setElement(A, val, 0, 2)
    LAGr_Matrix_setElement(A, val, 0, 3)
    LAGr_Matrix_setElement(A, val, 1, 2)
    LAGr_Matrix_setElement(A, val, 1, 3)
    LAGr_Matrix_setElement(A, val, 2, 4)
    LAGr_Matrix_setElement(A, val, 3, 5)
    // lower triangle
    LAGr_Matrix_setElement(A, val, 2, 0)
    LAGr_Matrix_setElement(A, val, 3, 0)
    LAGr_Matrix_setElement(A, val, 2, 1)
    LAGr_Matrix_setElement(A, val, 3, 1)
    LAGr_Matrix_setElement(A, val, 4, 2)
    LAGr_Matrix_setElement(A, val, 5, 3)

    LAGr_Matrix_new(&frontier, GrB_UINT64, n, bit_matrix_ncols)
    LAGr_Matrix_new(&next, GrB_UINT64, n, bit_matrix_ncols)
    LAGr_Matrix_new(&seen, GrB_UINT64, n, bit_matrix_ncols)
    LAGr_Matrix_new(&not_seen, GrB_UINT64, n, bit_matrix_ncols)
    LAGr_Matrix_new(&PopCount, GrB_UINT64, n, bit_matrix_ncols)

    LAGr_Vector_new(&popcount, GrB_UINT64, n)
    LAGr_Vector_new(&ones, GrB_UINT64, n)
    LAGr_Vector_new(&level_v, GrB_UINT64, n)
    LAGr_Vector_new(&sp, GrB_UINT64, n)

    // initialize frontier matrix
//    LAGr_Matrix_setElement(frontier, 1L << 63, 0, 0)
//    LAGr_Matrix_setElement(frontier, 1L << 62, 1, 0)
    // to compute closeness centrality, start off with the diagonal as a frontier
    create_diagonal_bit_matrix(frontier);

    // initialize ones vector
    create_ones_vector(ones);

        // initialize seen matrix with all explicit zeros...
    for (GrB_Index i = 0; i < n; i++) {
        for (GrB_Index j = 0; j < bit_matrix_ncols; j++) {
            LAGr_Matrix_setElement(seen, 0, i, j)
        }
    }
    // ...except where the traversal starts
    LAGr_eWiseAdd(seen, NULL, NULL, GrB_BOR_UINT64, seen, frontier, NULL)
    LAGr_apply(not_seen, NULL, NULL, GrB_BNOT_UINT64, seen, NULL)

    // traversal
    for (GrB_Index level = 1; level < n; level++) {
        printf("========================= Level %2ld =========================\n\n", level);
        // level_v += 1
        LAGr_eWiseAdd(level_v, NULL, NULL, GrB_PLUS_UINT64, level_v, ones, NULL)

        // next = A^T * frontier = A * frontier
        LAGr_mxm(next, NULL, NULL, GxB_BOR_BAND_UINT64, A, frontier, NULL)

        // next = next & ~seen // n.b. masking is not applicable
        LAGr_eWiseMult(next, NULL, NULL, GrB_BAND_UINT64, next, not_seen, NULL)

        GrB_Matrix_apply(PopCount, NULL, NULL, op_popcount, next, NULL);
        LAGr_reduce(popcount, NULL, NULL, GxB_PLUS_UINT64_MONOID, PopCount, NULL)
        LAGr_reduce(&total_popcount, NULL, GxB_PLUS_UINT64_MONOID, popcount, NULL)
        if (total_popcount == 0) {
            printf("no new vertices found\n");
            break;
        }

//        GrB_Matrix_apply(popcount, NULL, NULL, GxB_PLUS_UINT64_MONOID, PopCount, NULL);
        LAGr_eWiseMult(popcount, NULL, NULL, GrB_TIMES_UINT64, popcount, level_v, NULL)
        LAGr_eWiseAdd(sp, NULL, NULL, GrB_PLUS_UINT64, sp, popcount, NULL)

        // seen = seen | next
        LAGr_eWiseAdd(seen, NULL, NULL, GrB_BOR_UINT64, seen, next, NULL)
        LAGr_apply(not_seen, NULL, NULL, GrB_BNOT_UINT64, seen, NULL)

        print_bit_matrices(frontier, next, seen, popcount, sp);

        // frontier = next
        LAGr_Matrix_dup(&frontier, next)
    }

    LAGRAPH_FREE_ALL ;
    LAGraph_finalize ( ) ;
}

