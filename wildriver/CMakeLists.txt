if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  message("Configuring for MSVC")
  cmake_minimum_required(VERSION 3.4)
else()
  cmake_minimum_required(VERSION 2.8)
endif()

project(WILDRIVER)

if (DEFINED DEBUG)
  add_definitions(-DDEBUG=${DEBUG})
  set(CMAKE_BUILD_TYPE Debug)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
  message("Building a debug release")
else()
  set(CMAKE_BUILD_TYPE Release)
  add_definitions(-DNDEBUG=1)
  message("Building an optimized release")
endif()

if (STATIC)
  set(WILDRIVER_LIBRARY_TYPE STATIC)
else()
  set(WILDRIVER_LIBRARY_TYPE SHARED)
endif()

if (DEFINED DEVEL AND NOT DEVEL EQUAL 0)
  message("Development warnings enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")

  # g++ specific
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    message("Using g++ warnings.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmissing-declarations")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wstrict-aliasing")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")
    if ( (GCC_VERSION VERSION_GREATER 5.0 OR GCC_VERSION VERSION_EQUAL 5.0))
      message("Using g++ >= 5.0 warnings.")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
    endif()
  endif()
endif()


if (DEFINED WILDRIVER_DIMENSION_TYPE)
  add_definitions(-DWILDIRVER_DIMENSION_TYPE=${WILDRIVER_DIMENSION_TYPE})
endif()

if (DEFINED WILDRIVER_INDEX_TYPE)
  add_definitions(-DWILDIRVER_INDEX_TYPE=${WILDRIVER_INDEX_TYPE})
endif()

if (DEFINED WILDRIVER_VALUE_TYPE)
  add_definitions(-DWILDIRVER_VALUE_TYPE=${WILDRIVER_VALUE_TYPE})
endif()

# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
      message("WARNING: Attempting to substitute -std=c++0x for -std=c++11.")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Doverride=")
    else()
      # set c++11 for g++ >= 4.7
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()
else()
  # set c++11 standard for other compilers
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()


if (NOT WIN32)
  # use gnu directories
  include(GNUInstallDirs)
else()
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

if (MSVC_IDE)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
endif()

# must be before including src
if (DEFINED TESTS AND NOT TESTS EQUAL 0)
  enable_testing()
endif()

include_directories("include")
add_subdirectory("include")
add_subdirectory("src")


if (DEFINED EXAMPLES)
  add_subdirectory("examples")
endif()

