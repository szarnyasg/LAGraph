//------------------------------------------------------------------------------
// LAGraph_cdlp: community detection using label propagation
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

// LAGraph_cdlp: Contributed by Gabor Szarnyas and Balint Hegyi,
// Budapest University of Technology and Economics
// (with accented characters: G\'{a}bor Sz\'{a}rnyas and B\'{a}lint Hegyi,
// using LaTeX syntax). https://inf.mit.bme.hu/en/members/szarnyasg .

// ## Background
//
// This function was originally written for the LDBC Graphalytics benchmark.
//
// The community detection using label propagation (CDLP) algorithm is
// defined both for directed and undirected graphs.
//
// The definition implemented here is described in the following document:
// https://ldbc.github.io/ldbc_graphalytics_docs/graphalytics_spec.pdf
//
// The algorithm is based on the one given in the following paper:
//
// Usha Raghavan, Réka Albert, and Soundar Kumara. "Near linear time algorithm
// to detect community structures in large-scale networks". In: Physical
// Review E 76.3 (2007), p. 036106, https://arxiv.org/abs/0709.2938
//
// The key idea of the algorithm is that each vertex is assigned the label
// that is most frequent among its neighbors. To allow reproducible
// experiments, the algorithm is modified to guarantee deterministic behavior:
// it always picks the smallest label in case of a tie:
//
// min ( argmax_{l} (# of neighbors with label l) )
//
// In other words, we need to compute the *minimum mode value* (minmode) for
// the labels among the neighbors.
//
// For directed graphs, the definition is slightly more complicated:
// a neighbor that is connected through both an outgoing and on an incoming
// edge counts twice.
//
// ## Example (undirected)
//
// For an example, let's assume an undirected graph where vertex 1 has four
// neighbors {2, 3, 4, 5}, and the current labels in the graph are
// L = [3, 5, 4, 5, 4].
//
// In this example, the distribution of labels among the neighbors of vertex 1
// is {4 => 2, 5 => 2}, therefore, the minimum mode value is 4.
//
// Next, we capture this operation using GraphBLAS operations and
// data structures. Notice that the neighbors of vertex 1 are encoded
// as a sparse vector in the adjacency matrix:
//
// A = | 0 1 1 1 1 |
//     | 1 . . .   |
//     | 1 .       |
//     | 1 .       |
//     | 1         |
//
// To allow propagating the labels along edges, we use a diagonal matrix
// with the elements of the diagonal set to the values of L:
//
// diag(L) = | 3 0 0 0 0 |
//           | 0 5 0 0 0 |
//           | 0 0 4 0 0 |
//           | 0 0 0 5 0 |
//           | 0 0 0 0 4 |
//
// If we multiply the transposed of this matrix with diag(L), we get a matrix
// containing the labels of the neighbor nodes. (The transpose operation is
// required as we propagate the labels in the reverse direction of the edges).
//
// In the example, this gives:
//
// AL = A'*diag(L) = | 0 5 4 5 4 |
//                   | . . .     |
//
// Next, we needo compute the minimum mode value for each row. As it is
// difficult to capture this operation as a monoid, we extract each row of
// the matrix and compute the minimum mode outside of GraphBLAS terms.
//
// TODO: a possible optimization is that the first iteration is rather trivial
//
// * in the undirected case, each vertex gets the minimal initial label (=id)
//   of its neighbors.
// * in the directed case, each vertex gets the minimal initial label (=id)
//   of its neighbors which are doubly-linked (on an incoming and on an
//   outgoing edge). In the absence of such a neighbor, it picks the minimal
//   label of its neighbors (connected through either an incoming or through
//   an outgoing edge).

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&L) ;                 \
    if (sanitize) GrB_free (&S) ;   \
    GrB_free (&S) ;                 \
    GrB_free (&AL) ;                \
}

// TODO: remove this
void Print_Label_Matrix(GrB_Matrix m)
{
    GrB_Index row_vals;
    GrB_Matrix_nrows(&row_vals, m);

    printf("Label vec:\n");

    uint64_t value;

    printf(" ");
    for (GrB_Index i = 0; i < row_vals; i++)
    {
        printf(" %ld", i + 1);
    }
    printf("\n");

    printf("[");
    for (GrB_Index i = 0; i < row_vals; i++)
    {
        GrB_Matrix_extractElement(&value, m, i, i);
        printf(" %ld", value);
    }
    printf(" ]\n");
}

int cmp_uint64(const void *a, const void *b)  {
    if ( *(uint64_t *)a > *(uint64_t *)b ) return +1;
    if ( *(uint64_t *)a < *(uint64_t *)b ) return -1;
    return 0;
}

GrB_Info LAGraph_cdlp        // compute cdlp for all nodes in A
(
    GrB_Vector *CDLP_handle, // output vector
    const GrB_Matrix A,      // input matrix
    // TODO: handle the case when matrix is not symmetric
    bool symmetric,          // denote whether the matrix is symmetric
    bool sanitize,           // if true, ensure A is binary
    int itermax,             // max number of iterations,
    double *t                // t [0] = sanitize time, t [1] = cdlp time,
                             // in seconds
)
{
    GrB_Info info;

    // labels are stored as a diagonal matrix
    GrB_Matrix L;
    GrB_Matrix S;
    GrB_Matrix AL;
    GrB_Vector CDLP;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (CDLP_handle == NULL)
    {
        return GrB_NULL_POINTER;
    }

    // n = size of A (# of nodes in the graph)
    GrB_Index n;
    LAGRAPH_OK (GrB_Matrix_nrows(&n, A))
    //--------------------------------------------------------------------------
    // ensure input is binary and has no self-edges
    //--------------------------------------------------------------------------

    double tic [2] ;
    t [0] = 0 ;         // sanitize time
    t [1] = 0 ;         // CDLP time

    LAGraph_tic(tic);

    // TODO: maybe iteration can be terminated when L = L'
    if (sanitize)
    {
        LAGraph_tic (tic) ;

        // S = binary pattern of A
        LAGRAPH_OK (LAGraph_pattern (&S, A)) ;

        // remove all self edges
        LAGRAPH_OK (LAGraph_prune_diag (S)) ;
        t [0] = LAGraph_toc (tic) ;
    }
    else
    {
        // Use the input as-is, and assume it is binary with no self edges.
        // Results are undefined if this condition does not hold.
        S = A ;
    }

    LAGraph_tic (tic) ;

    GrB_Matrix_nrows(&n, A);

    // Initialize L with diagonal elements 1..n
    LAGRAPH_OK(GrB_Matrix_new(&L, GrB_UINT64, n, n))
    for (GrB_Index i = 0; i < n; i++)
    {
        LAGRAPH_OK(GrB_Matrix_setElement(L, i + 1, i, i))
    }

    GxB_Format_Value A_format = -1;
    LAGRAPH_OK (GxB_get (A, GxB_FORMAT, &A_format))
    if (A_format != GxB_BY_ROW)
    {
        LAGRAPH_ERROR(
            "CDLP algorithm only works on matrices stored by row (CSR)",
            GrB_INVALID_OBJECT
        )
    }

    for (int iteration = 0; iteration < itermax; iteration++)
    {
//        Print_Label_Matrix(L);
//        printf("\n\nIteration %d ----\n", iteration);

        // AL = A * L
        LAGRAPH_OK(GrB_Matrix_new(&AL, GrB_UINT64, n, n))
        LAGRAPH_OK(GrB_mxm(AL, NULL, NULL, GxB_PLUS_TIMES_UINT64, S, L, NULL ))
        GrB_free(&L);

        GrB_Index nz; // nz = # of nnz elements in the graph
        GrB_Index *I = NULL, *J = NULL; // for col/row indices of entries from AT
        uint64_t *V = NULL;

        LAGRAPH_OK (GrB_Matrix_nvals (&nz, AL));

        I = LAGraph_malloc (nz, sizeof(GrB_Index)) ;
        J = LAGraph_malloc (nz, sizeof(GrB_Index)) ;
        V = LAGraph_malloc (nz, sizeof(GrB_UINT64)) ;

        // initialize new matrix for L
        LAGRAPH_OK(GrB_Matrix_new(&L, GrB_UINT64, n, n))

        //--------------------------------------------------------------------------
        // extract tuples from AL, compute new labels
        //--------------------------------------------------------------------------
        LAGRAPH_OK (GrB_Matrix_extractTuples(I, J, V, &nz, AL));
        GrB_free(&AL);

        GrB_Index current_row = 0;
        GrB_Index row_start_index = 0;

        // We are iterating through the rows of the matrix, exploiting the
        // fact that SuiteSparse:GraphBLAS uses CSR by default and returns
        // the tuples in sorted order. Spec v3.0.1, GrB_Matrix_extractTuples:
        //
        // "The GraphBLAS API states the tuples can be returned in any order.
        // SuiteSparse:GraphBLAS chooses to always return them in sorted order,
        // depending on whether the matrix is stored by row or by column."
        //
        // Note that the column position (J[k]) is not relevant here (it
        // represents where the label 'comes from' which does not matter).
        for (GrB_Index k = 0; k <= nz; k++)
        {
            if (k == nz || I[k] > current_row)
            {
                // We find the minmode value using a sorted list.
                qsort(
                    &V[row_start_index],
                    k - row_start_index,
                    sizeof(uint64_t),
                    cmp_uint64
                );

                // Upon each change in value, we assess whether the current
                // run is longer than the previous longest run and if it is,
                // we change the longest run. This guarantees that we always
                // select the smallest one of the most frequent label.
                uint64_t max_value = V[row_start_index];
                uint64_t max_count = 1;
                uint64_t current_count = 1;

//                printf(" | %ld ", V[row_start_index]);
                for (int q = row_start_index + 1; q <= k; q++)
                {
//                    printf("%ld ", V[q]);
                    if (q < k && V[q] == V[q-1])
                    {
                        current_count++;
                    }
                    else // change of value
                    {
                        if (current_count > max_count)
                        {
                            max_value = V[q-1];
                            max_count = current_count;
                        }
                        current_count = 1;
                    }
                }
//                printf(" min_argmax=%ld \n", max_value);

                LAGRAPH_OK(GrB_Matrix_setElement(L, max_value, current_row, current_row))

                if (k < nz)
                {
                    row_start_index = k;
                    current_row = I[k];
                } // else the loop ends
            }
        }
        free(I);
        free(J);
        free(V);
    }
//    Print_Label_Matrix(L);

    //--------------------------------------------------------------------------
    // extract final labels to the result vector
    //--------------------------------------------------------------------------

    LAGRAPH_OK (GrB_Vector_new (&CDLP, GrB_UINT64, n))
    for (GrB_Index i = 0; i < n; i++)
    {
        uint64_t x;
        LAGRAPH_OK(GrB_Matrix_extractElement(&x, L, i, i))
        LAGRAPH_OK(GrB_Vector_setElement(CDLP, x, i))
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*CDLP_handle) = CDLP;
    CDLP = NULL;            // set to NULL so LAGRAPH_FREE_ALL doesn't free it
    LAGRAPH_FREE_ALL

    t[1] = LAGraph_toc(tic);
    return (GrB_SUCCESS);
}
