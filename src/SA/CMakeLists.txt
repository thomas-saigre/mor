# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

feelpp_add_application( sensitivity_analysis SRCS sensitivity_analysis.cpp
    PROJECT mor
    LINK_LIBRARIES OT Feelpp::feelpp_mor tbb # omp
)
