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

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&id2index) ;          \
}

#define ASSERT_TRUE(expr)                                   \
{                                                           \
    if(!(expr)) {                                           \
        fprintf(stderr, "Test failed: %s\nFile: %s:%d\n",   \
                #expr, __FILE__, __LINE__);                 \
        LAGRAPH_FREE_ALL;                                   \
        exit(EXIT_FAILURE);                                 \
    }                                                       \
}

uint64_t count_lines(const char* filepath) {
    uint64_t k = 0;
    char c;
    FILE* f = fopen(filepath, "r");
    while ((c = fgetc(f)) != EOF) {
        if (c == '\n') {
            k++;
        }
    }
    fclose(f);
    return k;
}

int main(int argc, char** argv) {

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Vector id2index = NULL;

    LAGRAPH_TRY_CATCH (LAGraph_init());
    
    bool remap = (strcmp (argv [1], "remap") == 0);

    uint64_t n = count_lines(argv[2]);
    uint64_t m = count_lines(argv[3]);

    FILE* fv = fopen(argv[2], "r");
    FILE* fe = fopen(argv[3], "r");

    GrB_Index* nodes = malloc(n * sizeof(GrB_Index));
    GrB_Index* srcs = malloc(m * sizeof(GrB_Index));
    GrB_Index* trgs = malloc(m * sizeof(GrB_Index));
    uint64_t * X = malloc(m * sizeof(uint64_t));

    // printf("vertices\n");
    // printf("================\n");

    GrB_Index node;
    GrB_Index i;

    i = 0;
    GrB_Index max = 0;
    while (fscanf (fv, "%ld\n", &node) != EOF)
    {
        // printf("%ld\n", node);
        nodes[i] = node;
        if (node > max) {
            max = node;
        }
        i++;
    }

    // printf("================\n");
    
    printf("nodes: %ld\n", n);
    printf("edges: %ld\n", m);
    printf("max: %ld\n", max);
    // printf("edges\n");
    // printf("================\n");

    GrB_Index src, trg;
    GrB_Index src_remapped, trg_remapped;
    double w;

    double tic1[2];
    LAGraph_tic(tic1);

    GrB_Index id_dimension;
    if (remap) {
        LAGRAPH_TRY_CATCH(LAGraph_dense_relabel(NULL, NULL, &id2index, nodes, n, NULL));
    }

    i = 0;
    while (fscanf (fe, "%ld %ld %lf\n", &src, &trg, &w) != EOF)
    {
        // printf("%ld: %d, %d\n", i, src, trg);
        if (remap) {
            LAGr_Vector_extractElement(&src_remapped, id2index, src);
            LAGr_Vector_extractElement(&trg_remapped, id2index, trg);
            srcs[i] = src_remapped;
            trgs[i] = trg_remapped;
        } else {
            srcs[i] = src;
            trgs[i] = trg;
        }

        X[i] = 1;
        i++;
    }

    GrB_Matrix A = NULL;
    if (remap) {
        GrB_Matrix_new(&A, GrB_UINT64, n, n);
    } else {
        GrB_Matrix_new(&A, GrB_UINT64, max+1, max+1);
    }
    GrB_Matrix_build(A, srcs, trgs, X, m, GrB_PLUS_UINT64);
    
    //GrB_Matrix_type(&type, )
    GxB_print(A, GxB_SUMMARY);

    double t1 = LAGraph_toc(tic1);

    double tic2[2];
    LAGraph_tic(tic2);
    GrB_Index ntri;
    LAGraph_tricount(&ntri, 4, 0, NULL, A);
    double t2 = LAGraph_toc(tic2);

    printf("remapping: %s\n", remap ? "true" : "false");
    printf("ntri: %d\n", ntri);
    printf("loading edges:   %12.6g sec\n", t1);
    printf("processing time: %12.6g sec\n", t2);

    LAGRAPH_FREE_ALL;
    LAGraph_finalize();
}
