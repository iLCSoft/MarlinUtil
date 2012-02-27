########################################################
# cmake file for building MarlinUtil
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( MarlinUtil )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 5 )
SET( ${PROJECT_NAME}_VERSION_PATCH 1 )


### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.1 REQUIRED ) # minimum required Marlin version
FIND_PACKAGE( CED 1.4 REQUIRED )
FIND_PACKAGE( CLHEP REQUIRED )
FIND_PACKAGE( GSL REQUIRED )

FOREACH( pkg Marlin CED CLHEP GSL )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

# export DEPENDS variables to MarlinUtilConfig.cmake
SET( MarlinUtil_DEPENDS_INCLUDE_DIRS ${CED_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS} )
SET( MarlinUtil_DEPENDS_LIBRARY_DIRS ${CED_LIBRARY_DIRS} ${CLHEP_LIBRARY_DIRS} )
SET( MarlinUtil_DEPENDS_LIBRARIES ${CED_LIBRARIES} ${CLHEP_LIBRARIES} )




### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC )

    FIND_PACKAGE( Doxygen )

    IF( DOXYGEN_EXECUTABLE )
        ADD_SUBDIRECTORY( ./doc )
    ELSE()
        MESSAGE( SEND_ERROR "Could not find doxygen required to build documentation" )
        MESSAGE( "Please install doxygen or set INSTALL_DOC to OFF" )
        MESSAGE( "" )
    ENDIF()

ENDIF()


### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi" ) # FIXME -pedantic
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( BEFORE ./source/include )
INSTALL_DIRECTORY( ./source/include/ DESTINATION ./include/marlinutil FILES_MATCHING PATTERN "*.h" )


# get list of all source files
AUX_SOURCE_DIRECTORY( ./source/src library_sources )
AUX_SOURCE_DIRECTORY( ./source/src/ann library_sources )
AUX_SOURCE_DIRECTORY( ./source/src/mille library_sources )

SET_SOURCE_FILES_PROPERTIES( "./source/src/ann/kd_pr_search.cpp" PROPERTIES COMPILE_FLAGS "-fno-strict-aliasing" )

ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )


# helper macro for adding new tests
# build with 'make tests'
ADD_CUSTOM_TARGET( tests )
MACRO( ADD_MARLINUTIL_TEST _name )
    # MarlinUtil test
    IF( BUILD_TESTING )
        ADD_EXECUTABLE( ${_name} "./source/tests/${_name}.cc" )
    ELSE()
        ADD_EXECUTABLE( ${_name} EXCLUDE_FROM_ALL "./source/tests/${_name}.cc" )
    ENDIF()
    ADD_DEPENDENCIES( tests ${_name} )
    TARGET_LINK_LIBRARIES( ${_name} ${PROJECT_NAME} )
    #INSTALL( TARGETS ${_name} DESTINATION bin )
ENDMACRO()

ADD_MARLINUTIL_TEST( testmarlinutil )


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( MarlinUtilConfig.cmake MarlinUtilConfigVersion.cmake MarlinUtilLibDeps.cmake )

