//------------------------------------------------------------------------------
// LAGraph/Test/TriangleCount/ttest7.c: test program for LAGraph_tricount
// with method 7
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

// Contributed by Gabor Szarnyas

// Usage: ttest7 ...

#include "LAGraph.h"


int main (int argc, char **argv)
{

    int i;
    printf("%d ...", i);
//    long N = 5000;
//    long i;
//    float* v1 = malloc(N*sizeof(long));
//    float* v2 = malloc(N*sizeof(long));
//    float* p = malloc(N*sizeof(long));
//
//    for (i=0; i<N; i++) {
//        v1[i] = 1.0f;
//        v2[i] = 2.0f;
//    }
//
//    #pragma omp target map(to: v1, v2) map(from: p)x
//    #pragma omp parallel for
//    for (i=0; i<N; i++) {
//        float res = 0.0f;
//        for (int j = 0; j < 100; j++) {
//            res += v1[i] * v2[i] / 3.12f * v1[i] / 8.3f;
//        }
//        p[i] = res;
//    }
//
//    float max_val;
//    #pragma omp parallel for reduction(max:max_val)
//    for (int idx = 0; idx < N; idx++) {
//        max_val = max_val > p[idx] ? max_val : p[idx];
//    }
//
//    free(v1);
//    free(v2);
//    free(p);

//for (i=0; i<N; i++) {printf("%02f\n", p[i]);}


}

