//------------------------------------------------------------------------------
// LAGraph/Test/MatrixExtractKeepDimensions/matrix_extract_keep_dimensions_test.c:
// test program for LAGraph_Matrix_extract_keep_dimensions
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

// Usage: matrix_extract_keep_dimensions_test can be used with binary input graphs
//
// matrix_extract_keep_dimensions_test binarymatrixfile.grb

#include "LAGraph.h"

uint64_t extract(const GrB_Matrix A, const GrB_Index i, const GrB_Index j) {
    uint64_t x;
    GrB_Info info = GrB_Matrix_extractElement(&x, A, i, j);
    if (info == GrB_NO_VALUE) {
        return -1;
    }
    return x;
}

void print_matrix(const GrB_Matrix A) {
    GrB_Index nrows, ncols;
    GrB_Matrix_nrows(&nrows, A);
    GrB_Matrix_ncols(&ncols, A);
    for (GrB_Index i = 0; i < nrows; i++) {
        printf("%4ld:", i);
        for (GrB_Index j = 0; j < ncols; j++) {
            uint64_t val = extract(A, i, j);
            if (val == -1) {
                printf("     ");
            } else {
                printf(" %4ld", val);
            }
        }
        printf("\n");
    }
    printf("\n");
}


#define LAGRAPH_FREE_ALL                            \
{                                                   \
    GrB_free (&A) ;                                 \
    GrB_free (&C) ;                                 \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Matrix A = NULL, C = NULL ;

    LAGraph_init ( ) ;
    int nthreads_max = LAGraph_get_nthreads ( ) ;
    if (nthreads_max == 0) nthreads_max = 1 ;

    //--------------------------------------------------------------------------
    // get the input matrix
    //--------------------------------------------------------------------------

    FILE *out = stdout ;

    FILE *f ;
    if (argc < 2)
    {
        printf ("Usage: matrix_extract_keep_dimensions_test matrix_file.mtx\n") ;
        return (GrB_INVALID_VALUE) ;
    }
    else
    {
        f = fopen (argv[1], "r") ;
        if (f == NULL)
        {
            printf ("unable to open file [%s]\n", argv[1]) ;
            return (GrB_INVALID_VALUE) ;
        }
    }
    LAGRAPH_OK (LAGraph_mmread (&A, f)) ;

    GrB_Index n;
    GrB_Matrix_nrows(&n, A);

    GrB_Matrix S = NULL;
    LAGraph_pattern(&S, A, GrB_NULL);

//    GxB_print(S, GxB_COMPLETE);
    print_matrix(S);
    LAGraph_reorder_vertices(&C, NULL, S, false);
    print_matrix(C);
//    GxB_print(C, GxB_COMPLETE);

    LAGRAPH_FREE_ALL ;
    LAGraph_finalize ( ) ;
}
