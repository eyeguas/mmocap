# ----------------------------------------------------------------------------
#   Basic Configuration
# ----------------------------------------------------------------------------
PROJECT(gummocap)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

set(PROJECT_VERSION "1.0.0")
string(REGEX MATCHALL "[0-9]" PROJECT_VERSION_PARTS "${PROJECT_VERSION}")
list(GET PROJECT_VERSION_PARTS 0 PROJECT_VERSION_MAJOR)
list(GET PROJECT_VERSION_PARTS 1 PROJECT_VERSION_MINOR)
list(GET PROJECT_VERSION_PARTS 2 PROJECT_VERSION_PATCH)
set(PROJECT_SOVERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")

#------------------------------------------------------
# Build type
#------------------------------------------------------

IF(NOT CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE "Debug" )
ENDIF()

#------------------------------------------------------
# Lib Names and Dirs
#------------------------------------------------------

if(WIN32)
    # Postfix of DLLs:
    set(PROJECT_DLLVERSION "${PROJECT_VERSION_MAJOR}${PROJECT_VERSION_MINOR}${PROJECT_VERSION_PATCH}")
    set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/../../winlibs")
    set(CMAKE_LIBRARY_PATH "${CMAKE_INSTALL_PREFIX}/lib")
    set(CMAKE_INCLUDE_PATH "${CMAKE_INSTALL_PREFIX}/include")

else()
    # Postfix of so's:
    set(PROJECT_DLLVERSION)
endif()
set(CMAKE_INCLUDE_DIRS_CONFIGCMAKE ${CMAKE_INSTALL_PREFIX}/include)
set(CMAKE_LIB_DIRS_CONFIGCMAKE ${CMAKE_INSTALL_PREFIX}/lib)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_INSTALL_PREFIX}/lib/cmake/)
message( STATUS "CMAKE_MODULE_PATH:             ${CMAKE_MODULE_PATH}")

# ----------------------------------------------------------------------------
#   Find Dependencies
# ----------------------------------------------------------------------------
FIND_PACKAGE(gurealea 	1.0.0	REQUIRED)
FIND_PACKAGE(gumocap 	1.0.0	REQUIRED)
FIND_PACKAGE(guscene 	1.0.0	REQUIRED)
FIND_PACKAGE(guai 	1.0.0	REQUIRED)
FIND_PACKAGE(gugl 	1.0.0	REQUIRED)
FIND_PACKAGE(OpenCV 		REQUIRED)
FIND_LIBRARY(gsl_LIBS       REQUIRED  NAMES gsl)
FIND_LIBRARY(gslcblas_LIBS REQUIRED  NAMES  gslcblas)
set(REQUIRED_LIBRARIES gomp  ${OpenCV_LIBS}   ${GL_LIBS} ${GLU_LIBS} ${glut_LIBS}  ${gsl_LIBS}  ${gslcblas_LIBS} ${gurealea_LIBS}  ${gumocap_LIBS} ${gugl_LIBS} ${guscene_LIBS}   ${guimage_LIBS} ${guai_LIBS}  ${gu_LIBS})
# ----------------------------------------------------------------------------
#   PROJECT CONFIGURATION
#   force some variables that could be defined in the command line to be written to cache
# ----------------------------------------------------------------------------
OPTION(BUILD_SHARED_LIBS 	"Set to ON to build static libraries" ON)

    OPTION(INSTALL_DOC 	"Set to OFF to skip build/install Documentation" OFF)
IF (WIN32)
    SET(RUNTIME_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Directory for dlls and binaries")
    SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Directory for binaries")
    SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Directory for dlls")
ENDIF ()

IF (INSTALL_DOC)
	INCLUDE("${PROJECT_SOURCE_DIR}/generateDoc.cmake")
	GENERATE_DOCUMENTATION(${PROJECT_SOURCE_DIR}/dox.in)
ENDIF()


# ----------------------------------------------------------------------------
#   Uninstall target, for "make uninstall"
# ----------------------------------------------------------------------------
CONFIGURE_FILE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

# ----------------------------------------------------------------------------
# create configuration file from .in file (If you use windows take care with paths)
# ----------------------------------------------------------------------------
if (NOT CMAKE_INSTALL_PREFIX)
    if (WIN32)
        set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR} CACHE INTERNAL "" FORCE)
    elseif()
        set(CMAKE_INSTALL_PREFIX "/usr" CACHE INTERNAL "" FORCE)
    endif()
endif()


set(CMAKE_INCLUDE_DIRS_CONFIGCMAKE ${CMAKE_INSTALL_PREFIX}/include/)
set(CMAKE_LIB_DIRS_CONFIGCMAKE ${CMAKE_INSTALL_PREFIX}/lib CACHE PATH "Output directory for libraries")

CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/config.cmake.in" "${PROJECT_BINARY_DIR}/Find${PROJECT_NAME}.cmake")
INSTALL(FILES "${PROJECT_BINARY_DIR}/Find${PROJECT_NAME}.cmake" DESTINATION lib/cmake/ )


IF (UNIX)
    CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/${PROJECT_NAME}.pc.in" "${PROJECT_BINARY_DIR}/${PROJECT_NAME}.pc" @ONLY INMEDIATE)
    INSTALL(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}.pc" DESTINATION lib/pkgconfig/ )
ENDIF()


# ----------------------------------------------------------------------------
#   Program Optimization and debug (Extracted from OpenCV)
# ----------------------------------------------------------------------------
set(WARNINGS_ARE_ERRORS 		OFF CACHE BOOL "Treat warnings as errors")
set(WHOLE_PROGRAM_OPTIMIZATION 	OFF CACHE BOOL "Flags for whole program optimization.")
# set (SUBREGION_XOR 1)
IF(SUBREGION_XOR)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DSUBREGION_XOR ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSUBREGION_XOR")

ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX OR MINGW)
    set(ENABLE_PROFILING 		OFF CACHE BOOL "Enable profiling in the GCC compiler (Add flags: -g -pg)")
    set(USE_OMIT_FRAME_POINTER 	ON CACHE BOOL "Enable -fomit-frame-pointer for GCC")
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES arm*) # We can use only -O2 because the -O3 causes gcc crash
        set(USE_O2 ON CACHE BOOL "Enable -O2 for GCC")
        set(USE_FAST_MATH ON CACHE BOOL "Enable -ffast-math for GCC")
    endif()
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES powerpc*)
        set(USE_O3 ON CACHE BOOL "Enable -O3 for GCC")
        set(USE_POWERPC ON CACHE BOOL "Enable PowerPC for GCC")
    endif ()
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES amd64* OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64*)
        set(USE_O3 ON CACHE BOOL "Enable -O3 for GCC")
        set(USE_FAST_MATH ON CACHE BOOL "Enable -ffast-math for GCC")
        set(USE_MMX ON CACHE BOOL "Enable MMX for GCC")
        set(USE_SSE ON CACHE BOOL "Enable SSE for GCC")
        set(USE_SSE2 ON CACHE BOOL "Enable SSE2 for GCC")
        set(USE_SSE3 ON CACHE BOOL "Enable SSE3 for GCC")
    endif()
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES i686* OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86)
        set(USE_O3 ON CACHE BOOL "Enable -O3 for GCC")
        set(USE_FAST_MATH ON CACHE BOOL "Enable -ffast-math for GCC")
        set(USE_MMX ON CACHE BOOL "Enable MMX for GCC")
        set(USE_SSE ON CACHE BOOL "Enable SSE for GCC")
        set(USE_SSE2 ON CACHE BOOL "Enable SSE2 for GCC")
        set(USE_SSE3 ON CACHE BOOL "Enable SSE3 for GCC")
    endif ()

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

    if(WARNINGS_ARE_ERRORS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
    endif()

    # The -Wno-long-long is required in 64bit systems when including sytem headers.
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64* OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES amd64*)
		  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long")
    endif()

    # Whole program optimization
    if(WHOLE_PROGRAM_OPTIMIZATION)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fwhole-program --combine")
    endif()

    # Other optimizations
    if(USE_OMIT_FRAME_POINTER)
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fomit-frame-pointer")
    endif()
    if(USE_O2)
       set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
    endif()
    if(USE_O3)
       set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    endif()
    if(USE_FAST_MATH)
       set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")
    endif()
    if(USE_POWERPC)
       set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mcpu=G3 -mtune=G5")
    endif()
    if(USE_MMX)
       set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mmmx")
    endif()
    if(USE_SSE)
       set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse")
    endif()
    if(USE_SSE2)
       set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse2")
    endif()
#     if(USE_SSE3 AND NOT MINGW) # SSE3 should be disabled under MingW because it generates compiler errors
    if(USE_SSE3 ) # changed by RMS
       set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse3")
    endif()

    if(ENABLE_PROFILING)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pg -g")
    else()
        if(NOT APPLE)
            set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffunction-sections")
        endif()
    endif()
    
    #add openmp support  
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -fopenmp ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -O0 -DDEBUG -D_DEBUG -W -Wextra -Wno-return-type -fopenmp -msse3")
ENDIF	()


# ----------------------------------------------------------------------------
#   PACKAGE GENERATOR
# ----------------------------------------------------------------------------

IF (UNIX)
SET(CPACK_GENERATOR "DEB")
SET(CPACK_PACKAGE_CONTACT "rmsalinas@uco.es")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "rmsalinas@uco.es")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Gu Markerless Mocap Library")
SET(CPACK_PACKAGE_VENDOR "Rafael Munoz-Salinas. AVA Group. University of Cordoba. Spain")
SET(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "gu (>= 2.0) ,  guimage (>=1.0.4), gurealea (>=1.0), gumocap (>=1.0), guscene(>=1.0), gustereo3 (>=1.0), guavi (>=1.2), guai (>=1.0), gugl (>=1.0) ")
SET(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
INCLUDE(CPack)
ENDIF()
# ----------------------------------------------------------------------------
#   PROJECT CONTENT
# ----------------------------------------------------------------------------
ADD_SUBDIRECTORY(src)
# ADD_SUBDIRECTORY(src/algorithms)
ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(tests)
# ----------------------------------------------------------------------------
# display status message for important variables
# ----------------------------------------------------------------------------
message( STATUS )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
message( STATUS "General configuration for ${PROJECT_NAME} ${VERSION}")
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
message( STATUS )
message( STATUS "Compiler:                  ${CMAKE_COMPILER}")
message( STATUS "C++ flags (Release):       ${CMAKE_CXX_FLAGS_RELEASE}")
message( STATUS "C++ flags (Debug):         ${CMAKE_CXX_FLAGS_DEBUG}")
message( STATUS "CMAKE_CXX_FLAGS=           ${CMAKE_CXX_FLAGS}")
message( STATUS "CMAKE_C_FLAGS=             ${CMAKE_C_FLAGS}")
message( STATUS "Linker flags (Release):    ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
message( STATUS "Linker flags (Debug):      ${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")

MESSAGE( STATUS )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR = ${CMAKE_SYSTEM_PROCESSOR}" )
MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}" )
MESSAGE( STATUS )
MESSAGE( STATUS "---------------------------     Documentation     -----------------------------" )
MESSAGE( STATUS )
MESSAGE( STATUS "INSTALL_DOC = ${INSTALL_DOC}" )
MESSAGE( STATUS "USE_LATEX = ${USE_LATEX}" )
MESSAGE( STATUS "USE_DOT = ${USE_DOT}" )
MESSAGE( STATUS "USE_CHM = ${USE_CHM}" )
MESSAGE( STATUS )
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS )

