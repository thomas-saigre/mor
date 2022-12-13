# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

feelpp_add_application( deterministic_sensitivity_analysis SRCS deterministic_sensitivity_analysis.cpp
    PROJECT mor
    LINK_LIBRARIES Feelpp::feelpp_mor tbb # omp
)
