# helper macro for adding new tests
# build with 'make tests'
ADD_CUSTOM_TARGET( tests )
MACRO( ADD_MARLINUTIL_TEST _name )
    # MarlinUtil test
    ADD_EXECUTABLE( ${_name} "${_name}.cc" )
    ADD_TEST( NAME ${_name} COMMAND ${_name})
    TARGET_LINK_LIBRARIES( ${_name} ${PROJECT_NAME} )
ENDMACRO()

ADD_MARLINUTIL_TEST( testmarlinutil )
