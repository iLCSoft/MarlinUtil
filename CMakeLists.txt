########################################################
# cmake file for building MarlinUtil
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( MarlinUtil )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 18 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )


### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.1 REQUIRED ) # minimum required Marlin version
FIND_PACKAGE( CED 1.4 REQUIRED )
FIND_PACKAGE( GSL 2.1 REQUIRED )

#set(PACKAGE_VERSION_COMPATIBLE True)
FIND_PACKAGE( CLHEP REQUIRED )

FOREACH( pkg Marlin CED CLHEP GSL )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( SYSTEM ${${pkg}_INCLUDE_DIRS} )
        #LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

#added by Thorben Quast
FIND_PACKAGE( DD4hep REQUIRED COMPONENTS DDRec)
INCLUDE_DIRECTORIES(BEFORE SYSTEM ${DD4hep_INCLUDE_DIRS} )
#LINK_LIBRARIES( ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES})
ADD_DEFINITIONS( ${DD4hep_DEFINITIONS} )
FIND_PACKAGE( ROOT REQUIRED )
INCLUDE_DIRECTORIES( SYSTEM  ${ROOT_INCLUDE_DIRS} )
#LINK_LIBRARIES( ${ROOT_LIBRARIES} )
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )


# export DEPENDS variables to MarlinUtilConfig.cmake
SET( MarlinUtil_DEPENDS_INCLUDE_DIRS ${CED_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS} ${DD4hep_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} )
SET( MarlinUtil_DEPENDS_LIBRARY_DIRS ${CED_LIBRARY_DIRS} ${CLHEP_LIBRARY_DIRS} ${DD4hep_LIBRARY_DIRS} ${ROOT_LIBRARY_DIRS} )

IF( TARGET CLHEP::CLHEP )
  get_property( XX_clhep_lib TARGET CLHEP::CLHEP PROPERTY LOCATION_${CMAKE_BUILD_TYPE} )
  #MESSAGE(INFO "************** export CLHEP_LIBRARIES: ${XX_clhep_lib}"  )
  SET( MarlinUtil_DEPENDS_LIBRARIES ${CED_LIBRARIES} ${XX_clhep_lib} ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES} ${ROOT_LIBRARIES} )
ELSE()
  SET( MarlinUtil_DEPENDS_LIBRARIES ${CED_LIBRARIES} ${CLHEP_LIBRARIES} ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES} ${ROOT_LIBRARIES} )
ENDIF()


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
#ADD_DEFINITIONS( "-Wall -ansi" ) # FIXME -pedantic
#ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( BEFORE ./source/include )
INSTALL_DIRECTORY( ./source/include/ DESTINATION ./include/marlinutil FILES_MATCHING PATTERN "*.h" PATTERN "*.ipp")

# get list of all source files
AUX_SOURCE_DIRECTORY( ./source/src library_sources )
AUX_SOURCE_DIRECTORY( ./source/src/mille library_sources )
AUX_SOURCE_DIRECTORY( ./source/ann ann_sources )

ADD_SHARED_LIBRARY( MarlinUtilAnn ${ann_sources} )
SET_TARGET_PROPERTIES( MarlinUtilAnn PROPERTIES COMPILE_FLAGS "-w" )
INSTALL_SHARED_LIBRARY( MarlinUtilAnn DESTINATION lib )

ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
  ${ROOT_LIBRARIES}
  ${DD4hep_COMPONENT_LIBRARIES}
  ${DD4hep_LIBRARIES}
  ${Marlin_LIBRARIES}
  ${CED_LIBRARIES}
  ${CLHEP_LIBRARIES}
  ${GSL_LIBRARIES}
  MarlinUtilAnn
  )

#AUX_SOURCE_DIRECTORY( ./source/src/ann ann_library_sources )
SET_SOURCE_FILES_PROPERTIES( "./source/src/ann/kd_pr_search.cpp" PROPERTIES COMPILE_FLAGS "-fno-strict-aliasing" )

#ADD_SHARED_LIBRARY( ${PROJECT_NAME}_ann ${ann_library_sources} )
#INSTALL_SHARED_LIBRARY( ${PROJECT_NAME}_ann DESTINATION lib )

option(USE_EXTERNAL_CATCH2 "Link against an external Catch2 v3 static library, otherwise build it locally" AUTO)
set(USE_EXTERNAL_CATCH2 AUTO CACHE STRING "Link against an external Catch2 v3 static library, otherwise build it locally")
set_property(CACHE USE_EXTERNAL_CATCH2 PROPERTY STRINGS AUTO ON OFF)

IF( BUILD_TESTING )
  ADD_SUBDIRECTORY(source/tests)
ENDIF()

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( MarlinUtilConfig.cmake MarlinUtilConfigVersion.cmake )
