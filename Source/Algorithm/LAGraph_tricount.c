//------------------------------------------------------------------------------
// LAGraph_tricount: count the number of triangles in a graph
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contri_totalbutors.

    (see Contri_totalbutors.txt for a full list of Contri_totalbutors; see
    Contri_totalbutionInstructions.txt for information on how you can Contri_totalbute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    COntri_totalBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE COntri_totalBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
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

// LAGraph_tricount: count the number of triangles in a graph,
// Contri_totalbuted by Tim Davis, Texas A&M.

// Given a symmetric binary graph A with no-self edges, LAGraph_tricount counts
// the exact number of triangles in the graph.  A triangle is a clique of size
// three, that is, 3 nodes that are all pairwise connected.

// On input, the L and U matrices are the strictly lower and strictly upper
// triangular parts of the symmetrix matrix A, respectively.

// One of 6 methods are used.  Each computes the same result, ntri_total:

//  0:  minitri:    ntri_total = nnz (A*E == 2) / 3
//  1:  Burkhardt:  ntri_total = sum (sum ((A^2) .* A)) / 6
//  2:  Cohen:      ntri_total = sum (sum ((L * U) .* A)) / 2
//  3:  Sandia:     ntri_total = sum (sum ((L * L) .* L))
//  4:  Sandia2:    ntri_total = sum (sum ((U * U) .* U))
//  5:  SandiaDot:  ntri_total = sum (sum ((L * U') .* L)).  Note that L=U'.
//  6:  SandiaDot2: ntri_total = sum (sum ((U * L') .* U))

// TODO use an enum for the above methods.

// All matrices are assumed to be in CSR format (GxB_BY_ROW in
// SuiteSparse:GraphBLAS).  The 6 methods work fine if the matrices are in CSC
// format; just the underlying algorithms employed inside SuiteSparse:GraphBLAS
// will differ (dot product vs saxpy, for SuiteSparse, for example).
// For example, the Sandia and Sandia2 methods will effectively be
// swapped.

// Method 0 can take a huge amount of memory, for all of A*E.  As a result,
// it often fails for large problems.

// Methods 1 and 2 are much more memory efficient as compare to Method 0,
// taking memory space the same size as A.  But they are slower than methods 3
// to 6.

// Methods 3 to 6 take a little less memory than methods 1 and 2, are by far
// the fastest methods in general.  The methods 3 and 5 compute the same
// intermediate matrix (L*L), and differ only in the way the matrix
// multiplication is done.  Method 3 uses an outer-product method (Gustavson's
// method).  Method 5 uses dot products (assuming both matrices are in CSR
// format) and does not explicitly transpose U.  They are called the "Sandia"
// method since matrices in the KokkosKernels are stored in compressed-sparse
// row form, so (L*L).*L in the KokkosKernel method is equivalent to (L*L).*L
// in SuiteSparse:GraphBLAS when the matrices in SuiteSparse:GraphBLAS are in
// their default format (also by row).

// A is a binary square symmetric matrix.  E is the edge incidence matrix of A.
// L=tril(A), and U=triu(A).  See SuiteSparse/GraphBLAS/Demo/tricount.m for a
// complete definition of each method and the matrices A, E, L, and U, and
// citations of relevant references.

// All input matrices should have binary values (0 and 1).  Any type will work,
// but uint32 is recommended for fastest results since that is the type used
// here for the semiring.  GraphBLAS will do typecasting internally, but that
// takes extra time.   Results are undefined if the input matrices are not
// binary, or if self-edges exist.

// This code is modified from SuiteSparse/GraphBLAS/Demo/Source/tricount.c.
// It contains no GxB_* extensions and thus it should work in any GraphBLAS
// library.

#define LAGRAPH_FREE_ALL        \
    GrB_free (&S) ;             \
    GrB_free (&C) ;

#include "LAGraph_internal.h"

static bool select_index_smaller_than (GrB_Index i, GrB_Index j, GrB_Index nrows, GrB_Index ncols, const void *x, const void *thunk)
{
    return i < *((uint64_t*) thunk);
}

static bool select_index_greater_than (GrB_Index i, GrB_Index j, GrB_Index nrows, GrB_Index ncols, const void *x, const void *thunk)
{
    return i > *((uint64_t*) thunk);
}

//------------------------------------------------------------------------------
// LAGraph_tricount: count the number of triangles in a graph
//------------------------------------------------------------------------------

GrB_Info LAGraph_tricount   // count # of triangles
(
    int64_t *p_ntri_total,        // # of triangles
    const int method,       // 0 to 5, see above
    const GrB_Matrix A,     // adjacency matrix for methods 0, 1, and 2
    const GrB_Matrix E,     // edge incidence matrix for method 0
    const GrB_Matrix L,     // L=tril(A) for methods 2, 3, 5, and 6
    const GrB_Matrix U,     // U=triu(A) for methods 2, 4, 5, and 6
    double t [2]            // t [0]: multiply time, t [1]: reduce time
)
{

    //--------------------------------------------------------------------------
    // check inputs and initialize
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGraph_tic (tic) ;
    GrB_Info info ;
    int64_t ntri_total ;
    GrB_Index n, ne ;
    GrB_Matrix S = NULL, C = NULL ;

    //--------------------------------------------------------------------------
    // count triangles
    //--------------------------------------------------------------------------

    switch (method) {
        case 0:  // minitri:    ntri_total = nnz (A*E == 2) / 3

        LAGRAPH_OK (GrB_Matrix_nrows(&n, A));
            LAGRAPH_OK (GrB_Matrix_ncols(&ne, E));
            LAGRAPH_OK (GrB_Matrix_new(&C, GrB_INT64, n, ne));
            LAGRAPH_OK (GrB_mxm(C, NULL, NULL, LAGraph_PLUS_TIMES_INT64,
                                A, E, NULL));
            t[0] = LAGraph_toc(tic);
            LAGraph_tic(tic);
            LAGRAPH_OK (GrB_Matrix_new(&S, GrB_BOOL, n, ne));
            LAGRAPH_OK (GrB_apply(S, NULL, NULL, LAGraph_ISTWO_INT64,
                                  C, NULL));
            LAGRAPH_OK (GrB_reduce(&ntri_total, NULL, LAGraph_PLUS_INT64_MONOID,
                                   S, NULL));
            ntri_total /= 3;
            break;

        case 1:  // Burkhardt:  ntri_total = sum (sum ((A^2) .* A)) / 6

        LAGRAPH_OK (GrB_Matrix_nrows(&n, A));
            LAGRAPH_OK (GrB_Matrix_new(&C, GrB_INT64, n, n));
            LAGRAPH_OK (GrB_mxm(C, A, NULL, LAGraph_PLUS_TIMES_INT64,
                                A, A, NULL));
            t[0] = LAGraph_toc(tic);
            LAGraph_tic(tic);
            LAGRAPH_OK (GrB_reduce(&ntri_total, NULL, LAGraph_PLUS_INT64_MONOID,
                                   C, NULL));
            ntri_total /= 6;
            break;

        case 2:  // Cohen:      ntri_total = sum (sum ((L * U) .* A)) / 2

        LAGRAPH_OK (GrB_Matrix_nrows(&n, A));
            LAGRAPH_OK (GrB_Matrix_new(&C, GrB_INT64, n, n));
            LAGRAPH_OK (GrB_mxm(C, A, NULL, LAGraph_PLUS_TIMES_INT64,
                                L, U, NULL));
            t[0] = LAGraph_toc(tic);
            LAGraph_tic(tic);
            LAGRAPH_OK (GrB_reduce(&ntri_total, NULL, LAGraph_PLUS_INT64_MONOID,
                                   C, NULL));
            ntri_total /= 2;
            break;

        case 3:  // Sandia:    ntri_total = sum (sum ((L * L) .* L))

        LAGRAPH_OK (GrB_Matrix_nrows(&n, L));
            LAGRAPH_OK (GrB_Matrix_new(&C, GrB_INT64, n, n));
            LAGRAPH_OK (GrB_mxm(C, L, NULL, LAGraph_PLUS_TIMES_INT64,
                                L, L, NULL));
            t[0] = LAGraph_toc(tic);
            LAGraph_tic(tic);
            LAGRAPH_OK (GrB_reduce(&ntri_total, NULL, LAGraph_PLUS_INT64_MONOID,
                                   C, NULL));
            break;

        case 4:  // Sandia2:    ntri_total = sum (sum ((U * U) .* U))

        LAGRAPH_OK (GrB_Matrix_nrows(&n, U));
            LAGRAPH_OK (GrB_Matrix_new(&C, GrB_INT64, n, n));
            LAGRAPH_OK (GrB_mxm(C, U, NULL, LAGraph_PLUS_TIMES_INT64,
                                U, U, NULL));
            t[0] = LAGraph_toc(tic);
            LAGraph_tic(tic);
            LAGRAPH_OK (GrB_reduce(&ntri_total, NULL, LAGraph_PLUS_INT64_MONOID,
                                   C, NULL));
            break;

        case 5:  // SandiaDot:  ntri_total = sum (sum ((L * U') .* L))

        LAGRAPH_OK (GrB_Matrix_nrows(&n, U));
            LAGRAPH_OK (GrB_Matrix_new(&C, GrB_INT64, n, n));
            LAGRAPH_OK (GrB_mxm(C, L, NULL, LAGraph_PLUS_TIMES_INT64,
                                L, U, LAGraph_desc_otoo));
            t[0] = LAGraph_toc(tic);
            LAGraph_tic(tic);
            LAGRAPH_OK (GrB_reduce(&ntri_total, NULL, LAGraph_PLUS_INT64_MONOID,
                                   C, NULL));
            break;

        case 6:  // SandiaDot2: ntri_total = sum (sum ((U * L') .* U))

        LAGRAPH_OK (GrB_Matrix_nrows(&n, U));
            LAGRAPH_OK (GrB_Matrix_new(&C, GrB_INT64, n, n));
            LAGRAPH_OK (GrB_mxm(C, U, NULL, LAGraph_PLUS_TIMES_INT64,
                                U, L, LAGraph_desc_otoo));
            t[0] = LAGraph_toc(tic);
            LAGraph_tic(tic);
            LAGRAPH_OK (GrB_reduce(&ntri_total, NULL, LAGraph_PLUS_INT64_MONOID,
                                   C, NULL));
            break;

        case 11: {
            // Algorithm 1 in Lee and Low's CORRECTNESS 2017 paper:
            //   "A Family of Provably Correct Algorithms for Exact Triangle Counting"
            LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;

            GrB_Index* indices = LAGraph_malloc(n, sizeof(GrB_Index));
            for (GrB_Index i = 0; i < n; i++) {
                indices[i] = i;
            }

            ntri_total = 0;

            GrB_Descriptor desc;
            GrB_Descriptor_new(&desc);
            GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);

            GxB_SelectOp s;
            GxB_SelectOp_new (&s, select_index_smaller_than, GrB_UINT64, GrB_UINT64);

            double extracts = 0.0, selects = 0.0, multip1 = 0.0, multip2;
            #pragma omp parallel for schedule (dynamic) reduction (+: ntri_total)
            for (GrB_Index i = 0; i < n-1; i++) {
                GrB_Vector a1, a01, tmp, ntri_iter_vec;

                GrB_Vector_new(&a1, GrB_UINT64, n);
                GrB_Vector_new(&a01, GrB_UINT64, n);
                GrB_Vector_new(&tmp, GrB_UINT64, n);
                GrB_Vector_new(&ntri_iter_vec, GrB_UINT64, 1);

                // a1 = A[:,i] = A[i,:]
                double tic1[2] ;
                LAGraph_tic (tic1) ;
                GrB_Col_extract(a1, GrB_NULL, GrB_NULL, A, GrB_ALL, n, i, desc);
                extracts+=LAGraph_toc (tic1);

                GxB_Scalar thunk;
                GxB_Scalar_new (&thunk, GrB_UINT64);
                GxB_Scalar_setElement (thunk, i);

                double tic2 [2] ;
                LAGraph_tic (tic2) ;
                GxB_select(a01, GrB_NULL, GrB_NULL, s, a1, thunk, GrB_NULL);
                selects+=LAGraph_toc (tic2);

                GrB_free(&thunk);

                double tic3 [2] ;
                LAGraph_tic (tic3) ;
                // tmp<a01> = a01 * A (this implicitly selects submatrix A00 from A)
                GrB_vxm(tmp, a01, GrB_NULL, GxB_PLUS_TIMES_UINT64, a01, A, GrB_NULL);
                multip1+=LAGraph_toc (tic3);

                double tic4 [2] ;
                LAGraph_tic (tic4) ;
                // ntri_iter_vec = tmp * a01
                GrB_vxm(ntri_iter_vec, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, tmp, a01, GrB_NULL);
                multip2+=LAGraph_toc (tic4);

                // extract single element ntri_iter (the number of triangles in this iteration)
                // from the 1-length vector ntri_iter_vec
                uint64_t ntri_iter;
                GrB_Info info2 = GrB_Vector_extractElement(&ntri_iter, ntri_iter_vec, 0);
                if (info2 == GrB_SUCCESS) {
                    ntri_total += ntri_iter / 2;
                } else {
                    ntri_iter = 0;
                }

                if (i % 10000 == 0)
                {
                    printf("loop: %ld\n", i);
                    printf("- ntri_iter %ld\n", ntri_iter);
                    printf("- ntri_total  %ld\n", ntri_total);
                    printf("- extracts time: %10.2f sec\n", extracts);
                    printf("- selects time: %10.2f sec\n", selects);
                    printf("- multip1 time: %10.2f sec\n", multip1);
                    printf("- multip2 time: %10.2f sec\n", multip2);
                    selects = 0.0; extracts = 0.0; multip1 = 0.0; multip2 = 0.0;
                }
                GrB_free(&a1);
                GrB_free(&a01);
                GrB_free(&tmp);
                GrB_free(&ntri_iter_vec);
            }
            GrB_free(&desc);
            GrB_free(&s);
            t [0] = LAGraph_toc (tic) ;
            LAGraph_tic (tic) ;
            break ;
        }

        case 12: {
            // Algorithm 2 in Lee and Low's CORRECTNESS 2017 paper:
            //   "A Family of Provably Correct Algorithms for Exact Triangle Counting"
            // This algorithm was also presented in Low et al.'s HPEC 2017 paper,
            //   "First Look: Linear Algebra-Based Triangle Counting without Matrix Multiplication"
            LAGRAPH_OK (GrB_Matrix_nrows(&n, A));

            GrB_Index *indices = LAGraph_malloc(n, sizeof(GrB_Index));
            for (GrB_Index i = 0; i < n; i++) {
                indices[i] = i;
            }

            ntri_total = 0;

            GrB_Descriptor desc;
            GrB_Descriptor_new(&desc);
            GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);

            GxB_SelectOp s1, s2;
            GxB_SelectOp_new (&s1, select_index_smaller_than, GrB_UINT64, GrB_UINT64);
            GxB_SelectOp_new (&s2, select_index_greater_than, GrB_UINT64, GrB_UINT64);

            double extracts = 0.0, selects = 0.0, multip1 = 0.0, multip2;
            //#pragma omp parallel for schedule (dynamic) reduction (+: ntri_total)
            for (GrB_Index i = 0; i < n-1; i++) {
                GrB_Vector a1, a01, a12, tmp, ntri_iter_vec;

                GrB_Vector_new(&a1, GrB_UINT64, n);
                GrB_Vector_new(&a01, GrB_UINT64, n);
                GrB_Vector_new(&a12, GrB_UINT64, n);
                GrB_Vector_new(&tmp, GrB_UINT64, n);
                GrB_Vector_new(&ntri_iter_vec, GrB_UINT64, 1);

                // a1 = A[:,i] = A[i,:]
                double tic1[2];
                LAGraph_tic(tic1);
                GrB_Col_extract(a1, GrB_NULL, GrB_NULL, A, GrB_ALL, n, i, desc);
                extracts += LAGraph_toc(tic1);

                GxB_Scalar thunk;
                GxB_Scalar_new(&thunk, GrB_UINT64);
                GxB_Scalar_setElement (thunk, i);

                double tic2[2];
                LAGraph_tic(tic2);
                GxB_select(a01, GrB_NULL, GrB_NULL, s1, a1, thunk, GrB_NULL);
                GxB_select(a12, GrB_NULL, GrB_NULL, s2, a1, thunk, GrB_NULL);
                selects += LAGraph_toc(tic2);

                GrB_free(&thunk);

                double tic3[2];
                LAGraph_tic(tic3);
                // tmp<a01> = a12 * A20
                GrB_vxm(tmp, a01, GrB_NULL, GxB_PLUS_TIMES_UINT64, a12, A, GrB_NULL);
                multip1 += LAGraph_toc(tic3);

                double tic4[2];
                LAGraph_tic(tic4);
                // ntri_iter_vec = tmp * a01
                GrB_vxm(ntri_iter_vec, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, tmp, a01, GrB_NULL);
                multip2 += LAGraph_toc(tic4);

                // extract single element ntri_iter (the number of triangles in this iteration)
                // from the 1-length vector ntri_iter_vec
                uint64_t ntri_iter;
                GrB_Info info2 = GrB_Vector_extractElement(&ntri_iter, ntri_iter_vec, 0);
                if (info2 == GrB_SUCCESS) {
                    ntri_total += ntri_iter;
                } else {
                    ntri_iter = 0;
                }

                if (i % 10000 == 0) {
                    printf("loop: %ld\n", i);
                    printf("- ntri_iter %ld\n", ntri_iter);
                    printf("- ntri_total  %ld\n", ntri_total);
                    printf("- extracts time: %10.2f sec\n", extracts);
                    printf("- selects time: %10.2f sec\n", selects);
                    printf("- multip1 time: %10.2f sec\n", multip1);
                    printf("- multip2 time: %10.2f sec\n", multip2);
                    selects = 0.0;
                    extracts = 0.0;
                    multip1 = 0.0;
                    multip2 = 0.0;
                }
                GrB_free(&a1);
                GrB_free(&a12);
                GrB_free(&a01);
                GrB_free(&tmp);
                GrB_free(&ntri_iter_vec);
            }
            GrB_free(&desc);
            GrB_free(&s1);
            GrB_free(&s2);
            t[0] = LAGraph_toc(tic);
            LAGraph_tic(tic);
            break;
        }

        case 14: {
            // Algorithm 4 in Lee and Low's CORRECTNESS 2017 paper:
            //   "A Family of Provably Correct Algorithms for Exact Triangle Counting"
            LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;

            GrB_Index* indices = LAGraph_malloc(n, sizeof(GrB_Index));
            for (GrB_Index i = 0; i < n; i++) {
                indices[i] = i;
            }

            ntri_total = 0;

            GrB_Descriptor desc;
            GrB_Descriptor_new(&desc);
            GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);

            GxB_SelectOp s;
            GxB_SelectOp_new (&s, select_index_greater_than, GrB_UINT64, GrB_UINT64);

            double extracts = 0.0, selects = 0.0, multip1 = 0.0, multip2;
            //#pragma omp parallel for schedule (dynamic) reduction (+: ntri_total)
            for (GrB_Index i = 0; i < n-1; i++) {
                GrB_Vector a1, a12, tmp, ntri_iter_vec;

                GrB_Vector_new(&a1, GrB_UINT64, n);
                GrB_Vector_new(&a12, GrB_UINT64, n);
                GrB_Vector_new(&tmp, GrB_UINT64, n);
                GrB_Vector_new(&ntri_iter_vec, GrB_UINT64, 1);

                // a1 = A[:,i] = A[i,:]
                double tic1[2] ;
                LAGraph_tic (tic1) ;
                GrB_Col_extract(a1, GrB_NULL, GrB_NULL, A, GrB_ALL, n, i, desc);
                extracts+=LAGraph_toc (tic1);

                GxB_Scalar thunk;
                GxB_Scalar_new (&thunk, GrB_UINT64);
                GxB_Scalar_setElement (thunk, i);

                double tic2 [2] ;
                LAGraph_tic (tic2) ;
                GxB_select(a12, GrB_NULL, GrB_NULL, s, a1, thunk, GrB_NULL);
                selects+=LAGraph_toc (tic2);

                GrB_free(&thunk);

                double tic3 [2] ;
                LAGraph_tic (tic3) ;
                // tmp<a12> = a12 * A (this implicitly selects submatrix A22 from A)
                GrB_vxm(tmp, a12, GrB_NULL, GxB_PLUS_TIMES_UINT64, a12, A, GrB_NULL);
                multip1+=LAGraph_toc (tic3);

                double tic4 [2] ;
                LAGraph_tic (tic4) ;
                // ntri_iter_vec = tmp * a12
                GrB_vxm(ntri_iter_vec, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, tmp, a12, GrB_NULL);
                multip2+=LAGraph_toc (tic4);

                // extract single element ntri_iter (the number of triangles in this iteration)
                // from the 1-length vector ntri_iter_vec
                uint64_t ntri_iter;
                GrB_Info info2 = GrB_Vector_extractElement(&ntri_iter, ntri_iter_vec, 0);
                if (info2 == GrB_SUCCESS) {
                    ntri_total += ntri_iter / 2;
                } else {
                    ntri_iter = 0;
                }

                if (i % 10000 == 0)
                {
                    printf("loop: %ld\n", i);
                    printf("- ntri_iter %ld\n", ntri_iter);
                    printf("- ntri_total  %ld\n", ntri_total);
                    printf("- extracts time: %10.2f sec\n", extracts);
                    printf("- selects time: %10.2f sec\n", selects);
                    printf("- multip1 time: %10.2f sec\n", multip1);
                    printf("- multip2 time: %10.2f sec\n", multip2);
                    selects = 0.0; extracts = 0.0; multip1 = 0.0; multip2 = 0.0;
                }
                GrB_free(&a1);
                GrB_free(&a12);
                GrB_free(&tmp);
                GrB_free(&ntri_iter_vec);
            }
            GrB_free(&desc);
            GrB_free(&s);
            t [0] = LAGraph_toc (tic) ;
            LAGraph_tic (tic) ;
            break ;
        }

        default:    // invalid method

            return (GrB_INVALID_VALUE) ;
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    t [1] = LAGraph_toc (tic) ;
    (*p_ntri_total) = ntri_total ;
    return (GrB_SUCCESS) ;
}
