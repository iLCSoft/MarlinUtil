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

#--- Unit test setup. Use existing Catch2 or fetch a suitable version and build it
if(CMAKE_CXX_STANDARD GREATER_EQUAL 20)
  set(CATCH2_MIN_VERSION 3.4)
else()
  set(CATCH2_MIN_VERSION 3.1)
endif()
if(USE_EXTERNAL_CATCH2)
  if (USE_EXTERNAL_CATCH2 STREQUAL AUTO)
    find_package(Catch2 ${CATCH2_MIN_VERSION})
  else()
    find_package(Catch2 ${CATCH2_MIN_VERSION} REQUIRED)
  endif()
endif()

if(NOT Catch2_FOUND)
  MESSAGE(STATUS "Fetching local copy of Catch2 library for unit-tests...")
  # Build Catch2 with the default flags, to avoid generating warnings when we
  # build it
  SET(CXX_FLAGS_CMAKE_USED ${CMAKE_CXX_FLAGS})
  SET(CMAKE_CXX_FLAGS ${CXX_FLAGS_CMAKE_DEFAULTS})
  INCLUDE(FetchContent)
  fETCHCONTENT_DECLARE(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.0.1
    )
  FETCHCONTENT_MAKEAVAILABLE(Catch2)
  set(CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras ${CMAKE_MODULE_PATH})

  # Disable clang-tidy on external contents
  SET_TARGET_PROPERTIES(Catch2 PROPERTIES CXX_CLANG_TIDY "")

  # Hack around the fact, that the include directories are not declared as
  # SYSTEM for the targets defined this way. Otherwise warnings can still occur
  # in Catch2 code when templates are evaluated (which happens quite a bit)
  GET_TARGET_PROPERTY(CATCH2_IF_INC_DIRS Catch2 INTERFACE_INCLUDE_DIRECTORIES)
  SET_TARGET_PROPERTIES(Catch2 PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CATCH2_IF_INC_DIRS}")

  # Reset the flags
  SET(CMAKE_CXX_FLAGS ${CXX_FLAGS_CMAKE_USED})
endif()

INCLUDE(Catch)

ADD_EXECUTABLE(unittests unittests/TestHelixClass.cpp)
TARGET_LINK_LIBRARIES(unittests PUBLIC ${PROJECT_NAME} PRIVATE Catch2::Catch2WithMain)
CATCH_DISCOVER_TESTS(unittests
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
