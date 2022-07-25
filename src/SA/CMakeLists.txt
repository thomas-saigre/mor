set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)
feelpp_add_application( sensitivity_analysis SRCS SA.cpp
    PROJECT mor
    LINK_LIBRARIES Feelpp::feelpp_toolbox_heat_lib OT
)
