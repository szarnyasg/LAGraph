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

    // prepare array of IDs
    srand(0);

    const GrB_Index nids = 1 * 1000 * 1000;
    GrB_Index* identifiers = malloc(nids * sizeof(GrB_Index));
    for (int i = 0; i < nids; i++) {
        identifiers[i] = rand() * UINT_MAX;
    }

    // build mappings
    GrB_Index id_dimension;
    double tic[2];
    LAGraph_tic (tic);
    LAGRAPH_TRY_CATCH(LAGraph_dense_relabel(&Id2index, &Index2id, &id2index, identifiers, nids, &id_dimension));
    double time = LAGraph_toc(tic);
    printf("Dense relabel time: %.2f\n", time);

    LAGRAPH_FREE_ALL;
    LAGraph_finalize();
}
