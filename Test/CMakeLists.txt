#-------------------------------------------------------------------------------
# LAGraph/Test/CMakeLists.txt:  cmake script for testing LAGraph
#-------------------------------------------------------------------------------

link_libraries(lagraph_static)

#-------------------------------------------------------------------------------
# Test for LCC
#-------------------------------------------------------------------------------

add_executable(LCC_test LCC/lcctest.c)
add_test(
        NAME LCC_dir_test
        COMMAND LCC_test ../../Data/ldbc-directed-example-unweighted.mtx
        WORKSPACE LCC
)
add_test(
        NAME LCC_undir_test
        COMMAND LCC_test ../../Data/ldbc-undirected-example-unweighted.mtx
        WORKSPACE LCC
)

#-------------------------------------------------------------------------------
# Test for CDLP
#-------------------------------------------------------------------------------

add_executable(CDLP_test CDLP/cdlptest.c)
add_test(
        NAME CDLP_common_dir_test
        COMMAND CDLP_test ../../Data/ldbc-directed-example-unweighted.mtx
        WORKSPACE CDLP
)
add_test(
        NAME CDLP_common_undir_test
        COMMAND CDLP_test ../../Data/ldbc-undirected-example-unweighted.mtx
        WORKSPACE CDLP
)
add_test(
        NAME CDLP_dir_test
        COMMAND CDLP_test ../../Data/ldbc-cdlp-directed-example.mtx
        WORKSPACE CDLP
)
add_test(
        NAME CDLP_undir_test
        COMMAND CDLP_test ../../Data/ldbc-cdlp-undirected-example.mtx
        WORKSPACE CDLP
)
