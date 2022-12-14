feelpp_add_application( test_gather_ot
    SRCS mpi_ot.cpp
    PROJECT mor
    LINK_LIBRARIES OT
)

feelpp_add_application( test_gather_vect
    SRCS mpi_gather.cpp
    PROJECT mor
    LINK_LIBRARIES OT
)

add_test (
    NAME test_gather_ot
    COMMAND mpirun -np 4 feelpp_mor_test_gather_ot
)

add_test (
    NAME test_gather_vect
    COMMAND mpirun -np 4 feelpp_mor_test_gather_vect
)