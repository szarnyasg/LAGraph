#-------------------------------------------------------------------------------
# LAGraph/Test/CMakeLists.txt:  cmake script for testing LAGraph
#-------------------------------------------------------------------------------

link_libraries(lagraph_static)

#-------------------------------------------------------------------------------
# Test for MSBFS
#-------------------------------------------------------------------------------

add_executable(MSBFS_test MSBFS/msbfstest.c)
add_test(
        NAME MSBFS_test
        COMMAND MSBFS_test
        WORKSPACE MSBFS
)
