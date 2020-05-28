//------------------------------------------------------------------------------
// LAGraph/Test/DenseRelabel/dense_relabel_test.c: test program for LAGraph_dense_relabel
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

// Contributed by Marton Elekes and Gabor Szarnyas, BME

// Usage:
//
// dense_relabel_test

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&Id2index) ;          \
    GrB_free (&Index2id) ;          \
    GrB_free (&id2index) ;          \
    GrB_free (&id_vec) ;            \
    GrB_free (&index_vec) ;         \
    GrB_free (&ref_id_vec) ;        \
    GrB_free (&ref_index_vec) ;     \
}

// https://www.geeksforgeeks.org/binary-search/
GrB_Index binarySearch(GrB_Index arr[], GrB_Index l, GrB_Index r, GrB_Index x)
{
    if (r >= l) {
        GrB_Index mid = l + (r - l) / 2;

        // If the element is present at the middle
        // itself
        if (arr[mid] == x)
            return mid;

        // If element is smaller than mid, then
        // it can only be present in left subarray
        if (arr[mid] > x)
            return binarySearch(arr, l, mid - 1, x);

        // Else the element can only be present
        // in right subarray
        return binarySearch(arr, mid + 1, r, x);
    }

    // We reach here when element is not
    // present in array
    return -1;
}


int main(void) {

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Matrix Id2index = NULL;
    GrB_Matrix Index2id = NULL;
    GrB_Vector id2index = NULL;
    GrB_Vector id_vec = NULL;
    GrB_Vector index_vec = NULL;
    GrB_Vector ref_id_vec = NULL;
    GrB_Vector ref_index_vec = NULL;

    LAGRAPH_TRY_CATCH (LAGraph_init());

    // TODO: set threads here
    uint64_t nthreads = 8;
    LAGraph_set_nthreads (nthreads) ;

    // prepare array of IDs
    srand(0);

//    const GrB_Index nnodes =  5;                  const GrB_Index nedges = 15;
    const GrB_Index nnodes =  3.6 * 1000 * 1000;   const GrB_Index nedges = 447 * 1000 * 1000;

    GrB_Index* vertex_ids = malloc(nnodes * sizeof(GrB_Index));
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (GrB_Index i = 0; i < nnodes; i++) {
        vertex_ids[i] = rand();
    }
    printf("\n");

    GrB_Index* edge_srcs = malloc(nedges * sizeof(GrB_Index));
    GrB_Index* edge_trgs = malloc(nedges * sizeof(GrB_Index));
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (GrB_Index j = 0; j < nedges; j++) {
        edge_srcs[j] = vertex_ids[rand() % nnodes];
        edge_trgs[j] = vertex_ids[rand() % nnodes];
    }

    double tic[2];

    // build id <-> index mapping
    LAGraph_tic (tic);
    LAGRAPH_TRY_CATCH(LAGraph_dense_relabel(NULL, NULL, &id2index, vertex_ids, nnodes, NULL));
    GrB_Index *I = LAGraph_malloc(nnodes, sizeof(GrB_Index));
    GrB_Index *X = LAGraph_malloc(nnodes, sizeof(GrB_Index));
    GrB_Index nnodes2;
    GrB_Vector_extractTuples(I, X, &nnodes2, id2index);
    double time1 = LAGraph_toc(tic);
    printf("Vertex relabel time: %.2f\n", time1);

    /*
    printf("I: ");
    for (GrB_Index i = 0; i < nnodes; i++) printf("%d, ", I[i]);
    printf("\n");

    printf("X: ");
    for (GrB_Index i = 0; i < nnodes; i++) printf("%d, ", X[i]);
    printf("\n");
     */

    // remap edges
    LAGraph_tic (tic);
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (GrB_Index j = 0; j < nedges; j++) {
        GrB_Index src_index, trg_index;
        //GrB_Vector_extractElement_UINT64(&src_index, id2index, edge_srcs[j]);
        //GrB_Vector_extractElement_UINT64(&trg_index, id2index, edge_trgs[j]);
        src_index = X[binarySearch(I, 0, nnodes, edge_srcs[j])];
        trg_index = X[binarySearch(I, 0, nnodes, edge_trgs[j])];
//        printf("%d -> %d ==> %d -> %d\n", edge_srcs[j], edge_trgs[j], src_index, trg_index);
    }
    double time2 = LAGraph_toc(tic);
    printf("Edge relabel time: %.2f\n", time2);

    LAGRAPH_FREE_ALL;
    LAGraph_finalize();
}
