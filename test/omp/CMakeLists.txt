# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# if(OpenMP_CXX_FOUND)
#     target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
# endif()


feelpp_add_application( test_omp
    SRCS test_omp.cpp
    PROJECT mor
    LINK_LIBRARIES OpenMP::OpenMP_CXX
)

add_test (
    NAME test_openmp
    COMMAND feelpp_mor_test_omp
)