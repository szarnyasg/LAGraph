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
// min ( argmax_{l} (#neighbors with label l) )
//
// In other words, we need to compute the *minimum mode value* (minmode) for
// the labels among the neighbors.
//
// For directed graphs, the definition is refined slightly: a label on a 
// neighbor that is connected through both an outgoing and on an incoming edge
// counts twice:
//
// min ( argmax_{l} (#incoming neighbors with l + #outgoing neighbors with l) )
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
// Next, we need to compute the minimum mode value for each row. As it is
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

typedef struct {
  GrB_Index number_of_rows;
  GrB_Index *restrict row_offset;
  GrB_Index *restrict row_length;
} RowMap;

//#define CDLP_OMP_FREE_ALL                                                      \
//{                                                                              \
//    GrB_free(&row);                                                            \
//}
//
//#define CDLP_OMP_ERROR(message, info)                                          \
//{                                                                              \
//    fprintf (stderr, "LAGraph error: %s\n[%d]\n%s\nFile: %s Line: %d\n",       \
//        message, info, GrB_error ( ), __FILE__, __LINE__) ;                    \
//    CDLP_OMP_FREE_ALL                                                          \
//}
//
//#define CDLP_OMP_CHECK(method)                                                 \
//{                                                                              \
//    for_info = method ;                                                        \
//    if (! (for_info == GrB_SUCCESS || info == GrB_NO_VALUE))                   \
//    {                                                                          \
//        CDLP_OMP_ERROR("", info);                                              \
//        continue;                                                              \
//    }                                                                          \
}

#define LAGRAPH_FREE_ALL                                                       \
{                                                                              \
    GrB_free (&L) ;                                                            \
    if (sanitize) GrB_free (&S) ;                                              \
    GrB_free (&S) ;                                                            \
    GrB_free (&SL) ;                                                           \
}

// TODO: remove this
void Print_Label_Matrix(GrB_Matrix m) {
    GrB_Index row_vals;
    GrB_Matrix_nrows(&row_vals, m);

    printf("Label vec:\n");

    uint64_t value;

    printf(" ");
    for (GrB_Index i = 0; i < row_vals; i++) {
        printf(" %ld", i + 1);
    }
    printf("\n");

    printf("[");
    for (GrB_Index i = 0; i < row_vals; i++) {
        GrB_Matrix_extractElement(&value, m, i, i);
        printf(" %ld", value);
    }
    printf(" ]\n");
}

GrB_Index Get_Min_Modus(
    const GrB_Index *labels,
    const GrB_Index num_vals
) {
    GrB_Index modus_label = labels[num_vals - 1];
    GrB_Index modus_occurrence = 1;
    GrB_Index current_occurrence = 1;

    for (int64_t i = num_vals - 2; i >= 0; i--) {
        if (labels[i] != labels[i + 1]) {
            if (current_occurrence >= modus_occurrence) {
                modus_label = labels[i + 1];
                modus_occurrence = current_occurrence;
            }
            current_occurrence = 1;
        } else {
            current_occurrence++;
        }
    }

    if (current_occurrence >= modus_occurrence) {
        modus_label = labels[0];
    }
    return modus_label;
}

void Build_Row_Map(
    const GrB_Index *row_index,
    const GrB_Index element_count,
    RowMap *row_map
) {
    GrB_Index start_index = 0;
    GrB_Index run_length = 1;
    GrB_Index current_index = *row_index;

    row_map->number_of_rows = 0;

    for (GrB_Index i = 1; i < element_count; i++) {
        if (current_index != row_index[i]) {
            row_map->row_offset[row_map->number_of_rows] = start_index;
            row_map->row_length[row_map->number_of_rows] = run_length;
            row_map->number_of_rows++;

            start_index = i;
            run_length = 1;
            current_index = row_index[i];
        } else {
            run_length++;
        }
    }

    row_map->row_offset[row_map->number_of_rows] = start_index;
    row_map->row_length[row_map->number_of_rows] = run_length;
    row_map->number_of_rows++;
}

int CMP_uint64(const void *a, const void *b) {
    if (*(uint64_t *) a > *(uint64_t *) b) return +1;
    if (*(uint64_t *) a < *(uint64_t *) b) return -1;
    return 0;
}

GrB_Info LAGraph_cdlp
    (
        GrB_Vector *CDLP_handle, // output vector
        const GrB_Matrix A,      // input matrix
        // TODO: handle the case when matrix is not symmetric
        bool symmetric,          // denote whether the matrix is symmetric
        bool sanitize,           // if true, ensure A is binary
        int itermax,             // max number of iterations,
        double *t                // t [0] = sanitize time, t [1] = cdlp time,
        // in seconds
    ) {
    GrB_Info info;

    // Diagonal label matrix
    GrB_Matrix L = NULL;
    // Source adjacency matrix
    GrB_Matrix S = NULL;
    // S*L matrix
    GrB_Matrix SL = NULL;
    // Result CDLP vector
    GrB_Vector CDLP = NULL;

    GrB_Descriptor replace_desc = NULL;

    // Arrays holding extracted tuples
    GrB_Index *restrict I = NULL;
    GrB_Index *restrict V = NULL;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (CDLP_handle == NULL) {
        return GrB_NULL_POINTER;
    }

    //--------------------------------------------------------------------------
    // ensure input is binary and has no self-edges
    //--------------------------------------------------------------------------

    double tic[2];
    t[0] = 0;         // sanitize time
    t[1] = 0;         // CDLP time

    // TODO: remove perf measurements
    double sub_tick[2];
    double init_time = 0.0;
    double mxm_time = 0.0;
    double alloc_time = 0.0;
    double extract_time = 0.0;
    double minmod_time = 0.0;
    double free_time = 0.0;

    LAGraph_tic(tic);

    // TODO: the iteration can be terminated when L[n] = L[n-1]
    if (sanitize) {
        LAGraph_tic(tic);

        // S = binary pattern of A
        LAGRAPH_OK (LAGraph_pattern(&S, A))
        // Remove all self edges
        LAGRAPH_OK (LAGraph_prune_diag(S))
        t[0] = LAGraph_toc(tic);
    } else {
        // Use the input as-is, and assume it is binary with no self edges.
        // Results are undefined if this condition does not hold.
        S = A;
    }

    LAGraph_tic(sub_tick);

    GxB_Format_Value A_format = -1;
    LAGRAPH_OK (GxB_get(A, GxB_FORMAT, &A_format))
    if (A_format != GxB_BY_ROW) {
        LAGRAPH_ERROR(
            "CDLP algorithm only works on matrices stored by row (CSR)",
            GrB_INVALID_OBJECT
        )
    }

    LAGRAPH_OK(GrB_Descriptor_new(&replace_desc))
    LAGRAPH_OK(GrB_Descriptor_set(replace_desc, GrB_OUTP, GrB_REPLACE))

    // n = size of A (# of nodes in the graph)
    GrB_Index n;
    LAGRAPH_OK (GrB_Matrix_nrows(&n, A))
    // nz = # of nnz elements in the graph
    GrB_Index nz;
    LAGRAPH_OK (GrB_Matrix_nvals(&nz, A))

    // Initialize L with diagonal elements 1..n
    LAGRAPH_OK(GrB_Matrix_new(&L, GrB_UINT64, n, n))
    for (GrB_Index i = 0; i < n; i++) {
        LAGRAPH_OK(GrB_Matrix_setElement(L, i + 1, i, i))
    }

    LAGRAPH_OK(GrB_Matrix_new(&SL, GrB_UINT64, n, n))
    init_time = LAGraph_toc(sub_tick);

    I = LAGraph_malloc(nz, sizeof(uint64_t));
    if (I == NULL) {
        LAGRAPH_ERROR("Cannot initialize I", GrB_NULL_POINTER)
    }

    V = LAGraph_malloc(nz, sizeof(uint64_t));
    if (V == NULL) {
        LAGRAPH_ERROR("Cannot initialize J", GrB_NULL_POINTER)
    }

    RowMap row_map = {
        0,
        LAGraph_malloc(n, sizeof(GrB_Index)),
        LAGraph_malloc(n, sizeof(GrB_Index)),
    };
    if (row_map.row_offset == NULL) {
        LAGRAPH_ERROR("Cannot initialize row_map", GrB_NULL_POINTER)
    }
    if (row_map.row_length == NULL) {
        LAGRAPH_ERROR("Cannot initialize run_length", GrB_NULL_POINTER)
    }

    for (int iteration = 0; iteration < itermax; iteration++) {

        LAGraph_tic(sub_tick);

        // SL = A * L
        LAGRAPH_OK(GrB_mxm(
            SL,
            GrB_NULL,
            GrB_NULL,
            GxB_PLUS_TIMES_UINT64,
            S,
            L,
            replace_desc
        ))

//        GxB_print(SL, GxB_COMPLETE);

        mxm_time += LAGraph_toc(sub_tick);
        LAGraph_tic(sub_tick);

        GrB_Index element_count = nz;
        LAGRAPH_OK(GrB_Matrix_extractTuples_UINT64(
            I,
            GrB_NULL,
            V,
            &element_count,
            SL
        ))

        if (iteration == 0) {
            Build_Row_Map(I, element_count, &row_map);
        }

        const int nthreads = LAGraph_get_nthreads();
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (
            GrB_Index row_index = 0;
            row_index < row_map.number_of_rows;
            row_index++) {

            GrB_Index row_offset = row_map.row_offset[row_index];
            GrB_Index row_length = row_map.row_length[row_index];

            GrB_Index column_index = I[row_offset];

            qsort(
                V + row_offset,
                row_length,
                sizeof(uint64_t),
                CMP_uint64
            );

            GrB_Index min_mod = Get_Min_Modus(
                V + row_offset,
                row_length
            );

//            printf("%lu - ", column_index);
//            for (GrB_Index i = 0; i < row_length; i++) {
//                printf("%lu ", V[row_offset + i]);
//            }
//            printf("- %lu\n", min_mod);

            GrB_Matrix_setElement(
                L,
                min_mod,
                column_index,
                column_index
            );
        }

//        Print_Label_Matrix(L);
    }

    //--------------------------------------------------------------------------
    // extract final labels to the result vector
    //--------------------------------------------------------------------------

    LAGraph_tic(sub_tick);

    LAGRAPH_OK (GrB_Vector_new(&CDLP, GrB_UINT64, n))
    for (GrB_Index i = 0; i < n; i++) {
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

    double cdlp_write = LAGraph_toc(sub_tick);

//    printf("\n");
//    printf("MXM time:     %14.6f sec\n", mxm_time);
//    printf("Alloc time:     %14.6f sec\n", alloc_time);
//    printf("Extract time:     %14.6f sec\n", extract_time);
//    printf("Minmod time:     %14.6f sec\n", minmod_time);
//    printf("Free time:     %14.6f sec\n", free_time);
//    printf("CDLP w time:     %14.6f sec\n", cdlp_write);
//    printf("\n");

    t[1] = LAGraph_toc(tic);
    return (GrB_SUCCESS);
}

