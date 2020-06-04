//------------------------------------------------------------------------------
// LAGraph/Test/extract/extract_test.c: test program for GrB*extractElement
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

// Contributed by Tim Davis, Texas A&M

// Usage:  extract_test

#include "LAGraph.h"

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

#define LAGRAPH_FREE_ALL    \
    GrB_free (&Xvec) ;

#define CHECK(ok)                                                   \
    if (!(ok))                                                      \
    {                                                               \
        fprintf (stderr, "fail: %s %d\n", __FILE__, __LINE__) ;     \
        printf ("fail: %s %d\n", __FILE__, __LINE__) ;              \
        abort ( ) ;                                                 \
    }

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Vector Xvec = NULL ;
    LAGraph_init ( ) ;         // start LAGraph

    char *library_date ;
    GxB_get (GxB_LIBRARY_DATE, &library_date) ;
    printf ("SuiteSparse:GraphBLAS %s\n", library_date) ;

    //--------------------------------------------------------------------------
    // construct a nearly-dense vector
    //--------------------------------------------------------------------------

    if (argc < 2)
    {
        printf("Usage: ./extract_test <mode>\n");
        printf("<mode>: e = extract, b = binary search\n");
        return -1;
    }

    bool use_extract_element = (argv[1][0] == 'e');

    double tic [2] ;
    LAGraph_tic (tic) ;

    GrB_Index n = 64 * 1024 * 1024 ;
    printf ("extract test: n = %lu\n", n) ;
    LAGr_Vector_new (&Xvec, GrB_UINT64, n) ;
    for (uint64_t k = 1 ; k < n ; k++)
    {
        LAGr_Vector_setElement (Xvec, k, k) ;
    }

    double t = LAGraph_toc (tic) ;
    printf ("set     time %12.6f n/sec %12.6f million\n", t, 1e-6 * ((double) n-1) / t) ;

    LAGraph_tic (tic) ;
    GrB_Index ignore ;
    LAGr_Vector_nvals (&ignore, Xvec) ;
    t = LAGraph_toc (tic) ;
    printf ("wait    time %12.6f n/sec %12.6f million\n", t, 1e-6 * ((double) n-1) / t) ;

    LAGraph_tic (tic) ;
    //GxB_print (Xvec, GxB_SHORT) ;
    t = LAGraph_toc (tic) ;
    printf ("check   time %12.6f n/sec %12.6f million\n", t, 1e-6 * ((double) n-1) / t) ;

    //--------------------------------------------------------------------------
    // test binary searches
    //--------------------------------------------------------------------------

    uint64_t x ;

    srand(0);
    GrB_Index* indices = LAGraph_malloc(n, sizeof(GrB_Index));
    for (uint64_t k = 1 ; k < n ; k++)
    {
        indices[k] = rand() % (n-1) + 1;
    }

    uint64_t sum = 0;
    if (use_extract_element)
    {
        for (uint64_t k = 1 ; k < n ; k++)
        {
            LAGr_Vector_extractElement (&x, Xvec, indices[k]) ;
            sum += x;
        }
    }
    else
    {
        GrB_Index n2 = n;
        GrB_Index* I = LAGraph_malloc(n, sizeof(GrB_Index));
        uint64_t* X = LAGraph_malloc(n, sizeof(uint64_t));
        LAGr_Vector_extractTuples(I, X, &n2, Xvec);
        for (uint64_t k = 1 ; k < n ; k++)
        {
            sum += X[binarySearch(I, 0, n, indices[k])];
        }
        free(I) ;
        free(X) ;
    }

    t = LAGraph_toc (tic) ;
    printf ("extract time %12.6f n/sec %12.6f million\n", t, 1e-6 * ((double) n-1) / t) ;

    printf("totally unusable value: %ld\n", sum);

    info = GrB_Vector_extractElement (&x, Xvec, 0) ;
    CHECK (info == GrB_NO_VALUE) ;

    free(indices) ;
    LAGRAPH_FREE_ALL ;
    LAGraph_finalize ( ) ;
    return (GrB_SUCCESS) ;
}

