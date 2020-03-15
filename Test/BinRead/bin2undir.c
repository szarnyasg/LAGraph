//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M University

// usage:
// mtx2bin infile.mtx outfile.grb

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&A) ;                 \
}

int main (int argc, char **argv)
{
    GrB_Info info ;
    GrB_Matrix A = NULL ;

    if (argc < 3)
    {
        LAGRAPH_ERROR ("Usage: bin2undir infile.grb outfile.grb",
            GrB_INVALID_VALUE) ;
    }

    printf ("infile:  %s\n", argv [1]) ;
    printf ("outfile: %s\n", argv [2]) ;

    LAGRAPH_OK (LAGraph_init ( )) ;
    LAGRAPH_OK (LAGraph_binread (&A, argv [1]) ) ;
    LAGRAPH_OK (GrB_eWiseAdd( A, NULL, NULL, GrB_LOR, A, A, GrB_DESC_T1 ))
    LAGRAPH_OK (LAGraph_binwrite (&A, argv [2], NULL) ) ;

    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS) ;
}
