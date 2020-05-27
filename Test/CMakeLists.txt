#-------------------------------------------------------------------------------
# LAGraph/Test/CMakeLists.txt:  cmake script for testing LAGraph
#-------------------------------------------------------------------------------

link_libraries(lagraph_static)

#-------------------------------------------------------------------------------
# Test for DenseRelabel
#-------------------------------------------------------------------------------

add_executable(dense_relabel_test DenseRelabel/dense_relabel_test.c)
add_test(
        NAME DRL_test
        COMMAND dense_relabel_test
        WORKSPACE DRL
)
