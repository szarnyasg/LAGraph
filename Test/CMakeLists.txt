#-------------------------------------------------------------------------------
# LAGraph/Test/CMakeLists.txt:  cmake script for testing LAGraph
#-------------------------------------------------------------------------------

link_libraries(lagraph_static)

#-------------------------------------------------------------------------------
# Test for Matrix_extract_keep_dimensions
#-------------------------------------------------------------------------------

add_executable(LAGraph_reorder_vertices MatrixExtractKeepDimensions/matrix_extract_keep_dimensions_test.c)
add_test(
        NAME LAGraph_reorder_vertices_test
        COMMAND LAGraph_reorder_vertices_test ../../Data/ldbc-undirected-example.mtx
        WORKSPACE LAGraph_reorder_vertices
)
