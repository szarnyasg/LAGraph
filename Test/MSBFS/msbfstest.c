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
        return 0xcccccccccccccccc;
    }
    return x;
}

uint64_t extract_v(const GrB_Vector v, const GrB_Index i) {
    uint64_t x;
    GrB_Info info = GrB_Vector_extractElement(&x, v, i);
    if (info == GrB_NO_VALUE) {
        return 9999999;
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

void print_bit_matrices(const GrB_Matrix frontier, const GrB_Matrix next, const GrB_Matrix seen, const GrB_Vector next_popcount, const GrB_Vector sp) {
    GrB_Index nrows, ncols;
    GrB_Matrix_nrows(&nrows, frontier);
    GrB_Matrix_ncols(&ncols, frontier);
    printf("          frontier             next               seen         next_popcount    sp\n");
    for (GrB_Index i = 0; i < nrows; i++) {
        printf("%4ld:", i);
        for (GrB_Index j = 0; j < ncols; j++) {
            printf(
              " %016lx   %016lx   %016lx   %13ld   %3ld",
              extract(frontier, i, j), extract(next, i, j), extract(seen, i, j), extract_v(next_popcount, i), extract_v(sp, i));
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
//    X = repeat {b100..., b010..., b001..., ..., b...001} until we have n elements
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

uint64_t next_popcount(uint64_t x)
{
    int c = 0;
    for (; x != 0; x >>= 1)
        if (x & 1)
            c++;
    return c;
}

void fun_sum_popcount (void *z, const void *x)
{
    (*((uint64_t *) z))  = next_popcount(* ((uint64_t *) x));
}


int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Matrix A = NULL, frontier = NULL, next = NULL, seen = NULL, not_seen = NULL, Seen_PopCount = NULL ;

    GrB_Matrix Next_PopCount;
    GrB_Vector next_popcount;

    GrB_Vector ones, n_minus_one, level_v, sp, compsize, ccv;

    LAGraph_init ( ) ;

    // initializing unary operator for next_popcount
    GrB_UnaryOp op_popcount = NULL ;
    LAGRAPH_OK (GrB_UnaryOp_new(&op_popcount, fun_sum_popcount, GrB_UINT64, GrB_UINT64))
    GrB_Semiring semiring_bor_first = NULL, semiring_bor_second = NULL ;
    LAGRAPH_OK (GrB_Semiring_new(&semiring_bor_first, GxB_BOR_UINT64_MONOID, GrB_FIRST_UINT64))
    LAGRAPH_OK (GrB_Semiring_new(&semiring_bor_second, GxB_BOR_UINT64_MONOID, GrB_SECOND_UINT64))

    // create the input matrix
    const GrB_Index n = 8;
    const GrB_Index bit_matrix_ncols = (n+63)/64;

    LAGr_Matrix_new(&A, GrB_BOOL, n, n)
    const bool val = true;

    // upper triangle
    LAGr_Matrix_setElement(A, val, 0, 2)
    LAGr_Matrix_setElement(A, val, 0, 3)
    LAGr_Matrix_setElement(A, val, 1, 2)
    LAGr_Matrix_setElement(A, val, 1, 3)
    LAGr_Matrix_setElement(A, val, 2, 4)
    LAGr_Matrix_setElement(A, val, 3, 5)
    LAGr_Matrix_setElement(A, val, 6, 7)
    // lower triangle
    GrB_eWiseAdd(A, NULL, NULL, GrB_PLUS_UINT64, A, A, LAGraph_desc_otoo);

    LAGr_Matrix_new(&frontier, GrB_UINT64, n, bit_matrix_ncols)
    LAGr_Matrix_new(&next, GrB_UINT64, n, bit_matrix_ncols)
    LAGr_Matrix_new(&seen, GrB_UINT64, n, bit_matrix_ncols)
    LAGr_Matrix_new(&not_seen, GrB_UINT64, n, bit_matrix_ncols)
    LAGr_Matrix_new(&Next_PopCount, GrB_UINT64, n, bit_matrix_ncols)
    LAGr_Matrix_new(&Seen_PopCount, GrB_UINT64, n, bit_matrix_ncols)

    LAGr_Vector_new(&next_popcount, GrB_UINT64, n)
    LAGr_Vector_new(&ones, GrB_UINT64, n)
    LAGr_Vector_new(&n_minus_one, GrB_UINT64, n)
    LAGr_Vector_new(&level_v, GrB_UINT64, n)
    LAGr_Vector_new(&sp, GrB_UINT64, n)
    LAGr_Vector_new(&compsize, GrB_UINT64, n)
    LAGr_Vector_new(&ccv, GrB_FP64, n)

    // initialize frontier and seen matrices: to compute closeness centrality, start off with a diagonal
    create_diagonal_bit_matrix(frontier);
    LAGr_Matrix_dup(&seen, frontier)

    // initialize vectors
    GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL);
    GrB_assign(n_minus_one, NULL, NULL, n-1, GrB_ALL, n, NULL);

    // initialize

    // traversal
    for (GrB_Index level = 1; level < n; level++) {
        printf("========================= Level %2ld =========================\n\n", level);
        // level_v += 1
        LAGr_eWiseAdd(level_v, NULL, NULL, GrB_PLUS_UINT64, level_v, ones, NULL)

        // next = A^T * frontier = A * frontier
        LAGr_mxm(next, NULL, NULL, semiring_bor_second, A, frontier, NULL)

        LAGr_apply(not_seen, NULL, NULL, GrB_BNOT_UINT64, seen, NULL)
        // next = next & ~seen
        // we need to use eWiseAdd to see the union of value
        // the code previously dropped zero elements
        // DO NOT do that as it will render seen[i] = 0000 (implicit value)
        // and seen[j] = 1111 equivalent in not_seen
        LAGr_eWiseAdd(next, NULL, NULL, GrB_BAND_UINT64, next, not_seen, NULL)


        LAGr_apply(Next_PopCount, NULL, NULL, op_popcount, next, NULL);
        LAGr_reduce(next_popcount, NULL, NULL, GxB_PLUS_UINT64_MONOID, Next_PopCount, NULL)

        GrB_Index next_nvals;
        LAGr_Vector_nvals(&next_nvals, next_popcount)
        if (next_nvals == 0) {
            printf("no new vertices found\n");
            break;
        }
        // seen = seen | next
        LAGr_eWiseAdd(seen, NULL, NULL, GrB_BOR_UINT64, seen, next, NULL)

        // sp += (next_popcount * level)
        //   next_popcount * level is expressed as next_popcount *= level_v
        LAGr_eWiseMult(next_popcount, NULL, NULL, GrB_TIMES_UINT64, next_popcount, level_v, NULL)
        LAGr_eWiseAdd(sp, NULL, NULL, GrB_PLUS_UINT64, sp, next_popcount, NULL)

        print_bit_matrices(frontier, next, seen, next_popcount, sp);

        // frontier = next
        LAGr_Matrix_dup(&frontier, next)
    }
    // compsize = reduce(seen, row -> popcount(row))
    GrB_Matrix_apply(Seen_PopCount, NULL, NULL, op_popcount, seen, NULL);
    LAGr_reduce(compsize, NULL, NULL, GxB_PLUS_UINT64_MONOID, Seen_PopCount, NULL)

    // compute the closeness centrality value:
    //
    //          (C(p)-1)^2
    // CCV(p) = ----------
    //          (n-1)*s(p)

    // C(p)-1
    GrB_eWiseAdd(compsize, NULL, NULL, GrB_MINUS_UINT64, compsize, ones, NULL);
    // (C(p)-1)^2
    GrB_eWiseMult(compsize, NULL, NULL, GrB_TIMES_UINT64, compsize, compsize, NULL);

    GrB_eWiseMult(sp, NULL, NULL, GrB_TIMES_UINT64, n_minus_one, sp, NULL);
    GrB_eWiseMult(ccv, NULL, NULL, GrB_DIV_FP64, compsize, sp, NULL);

//    GxB_print(compsize, GxB_SHORT);
//    GxB_print(sp, GxB_SHORT);
    GxB_print(ccv, GxB_SHORT);

    LAGRAPH_FREE_ALL ;
    LAGraph_finalize ( ) ;
}
