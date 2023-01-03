feelpp_add_application( test_ot
    SRCS test_openturns.cpp
    PROJECT mor
    LINK_LIBRARIES OT
)

feelpp_add_application( polynomial_chaos
    SRCS polynomial_chaos.cpp
    PROJECT mor
    LINK_LIBRARIES OT
)

set(TESTFILES test.sh sobol.py)
install(FILES ${TESTFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


add_test (
    NAME test_openturns
    COMMAND ${CMAKE_COMMAND} -E env
    # ${PYTHON_EXECUTABLE} -m pytest test_eim.py -s -vv -k ${test}
    bash test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)