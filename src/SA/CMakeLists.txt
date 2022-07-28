set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)
feelpp_add_application( sensitivity_analysis SRCS SA.cpp
    PROJECT mor
    LINK_LIBRARIES OT Feelpp::feelpp_mor omp
)
