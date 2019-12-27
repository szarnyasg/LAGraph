#-------------------------------------------------------------------------------
# LAGraph/Test/CMakeLists.txt:  cmake script for testing LAGraph
#-------------------------------------------------------------------------------

link_libraries(lagraph_static)

#-------------------------------------------------------------------------------
# Test for TTest
#-------------------------------------------------------------------------------

add_executable(ttest TriangleCount/ttest.c)
add_test(
        NAME ttest
        COMMAND ttest ../../Data/ldbc-undirected-example.mtx
        WORKSPACE TriangleCount
)

#-------------------------------------------------------------------------------
# Test for TTest7
#-------------------------------------------------------------------------------

add_executable(ttest7 TriangleCount/ttest7.c)
add_test(
        NAME ttest7
        COMMAND ttest7 ../../Data/ldbc-undirected-example.mtx
        WORKSPACE TriangleCount
)
