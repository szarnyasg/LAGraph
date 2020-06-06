#-------------------------------------------------------------------------------
# LAGraph/Test/CMakeLists.txt:  cmake script for testing LAGraph
#-------------------------------------------------------------------------------

link_libraries(lagraph_static)

#-------------------------------------------------------------------------------
# Tests
#-------------------------------------------------------------------------------

add_executable(MSBFS_test MSBFS/msbfstest.c)
add_test(
        NAME MSBFS_test
        COMMAND MSBFS_test
        WORKSPACE MSBFS
)

add_executable(BIDIRBFS_test MSBFS/bidirbfstest.c)
add_test(
        NAME BIDIRBFS_test
        COMMAND BIDIRBFS_test
        WORKSPACE BIDIRBFS
)
