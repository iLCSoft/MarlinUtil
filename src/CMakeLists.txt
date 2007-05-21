######################################
# cmake file for building MarlinUtil
# @author Jan Engels, DESY
######################################

####################################################################################################
### INCLUDES
####################################################################################################

# include directories
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/include" )

# install include files
INSTALL( DIRECTORY "${${PROJECT_NAME}_SOURCE_DIR}/include"
		DESTINATION .
		PATTERN "*~" EXCLUDE
		PATTERN "*CVS*" EXCLUDE )


####################################################################################################
### SOURCES
####################################################################################################

# require proper c++
# FIXME -pedantic flag crashes under sl4
ADD_DEFINITIONS( "-Wall -ansi" )

# binary sources: list of all sources from which you want to create a binary
#SET( bin_srcs testmarlinutil.cc )

# exclude this sources from the the library
SET( exc_srcs testmarlinutil.cc )

# get list of all source files in this directory
AUX_SOURCE_DIRECTORY( . all_srcs )
AUX_SOURCE_DIRECTORY( ann ann_srcs )

# remove "./" from the list of source files
STRING( REPLACE "./" "" all_srcs "${all_srcs}" )

# exclude sources: filter out the sources specified above from the list of all sources
IF( bin_srcs OR exc_bin_srcs )
	LIST( REMOVE_ITEM all_srcs ${bin_srcs} ${exc_srcs} )
ENDIF()

# library sources
SET( lib_srcs ${all_srcs} ${ann_srcs} )


####################################################################################################
### DEPENDENCIES
####################################################################################################

# load macro
IF( DEFINED ${PROJECT_NAME}_DEPENDS OR DEFINED BUILD_WITH OR DEFINED LINK_WITH )
    # load macro
    IF( NOT EXISTS "${CMAKE_MODULE_PATH}/MacroCheckDeps.cmake" )
        MESSAGE( FATAL_ERROR "\nSorry, could not find MacroCheckDeps.cmake...\n"
            "Please set CMAKE_MODULE_PATH correctly with: "
            "cmake -DCMAKE_MODULE_PATH=<path_to_cmake_modules>" )
    ENDIF()
    INCLUDE( "${CMAKE_MODULE_PATH}/MacroCheckDeps.cmake" )
    CHECK_DEPS()
ENDIF()


####################################################################################################
### LIBRARY
####################################################################################################

# create library
ADD_LIBRARY( lib_${PROJECT_NAME} ${lib_srcs} )

# create symbolic lib target for calling target lib_XXX
ADD_CUSTOM_TARGET( lib DEPENDS lib_${PROJECT_NAME} )

# change lib_target properties
SET_TARGET_PROPERTIES( lib_${PROJECT_NAME} PROPERTIES
	# create *nix style library versions + symbolic links
	VERSION ${${PROJECT_NAME}_VERSION}
	SOVERSION ${${PROJECT_NAME}_SOVERSION}
	# allow creating static and shared libs without conflicts
	CLEAN_DIRECT_OUTPUT 1
	# avoid conflicts between library and binary target names
	OUTPUT_NAME ${PROJECT_NAME}	)

# install library
INSTALL( TARGETS lib_${PROJECT_NAME} DESTINATION lib PERMISSIONS
		OWNER_READ OWNER_WRITE OWNER_EXECUTE
		GROUP_READ GROUP_EXECUTE
		WORLD_READ WORLD_EXECUTE )


####################################################################################################
### BINS
####################################################################################################

IF( bin_srcs )
	# create symbolic bin target for calling targets bin_XXX
	ADD_CUSTOM_TARGET( bin )

	# create bin targets (remove .cc extensions from binary sources)
	FOREACH( src ${bin_srcs} )
		STRING( REPLACE ".cc" "" bin ${src} )
		LIST( APPEND bins ${bin} )
    ENDFOREACH()

	# loop through list of binaries
	FOREACH( bin ${bins} )
		# create target for binary
		ADD_EXECUTABLE( bin_${bin} EXCLUDE_FROM_ALL ${bin}.cc ${ext_bin_srcs} )
		# tell bin target to call this bin_XXX target
		ADD_DEPENDENCIES( bin bin_${bin} )
		# avoid conflicts between library and binary target names
		SET_TARGET_PROPERTIES( bin_${bin} PROPERTIES OUTPUT_NAME ${bin} )
		# link against project library
		TARGET_LINK_LIBRARIES( bin_${bin} lib_${PROJECT_NAME} )
		# install binary commando
		INSTALL( TARGETS bin_${bin} DESTINATION bin )
	ENDFOREACH()
ENDIF()

