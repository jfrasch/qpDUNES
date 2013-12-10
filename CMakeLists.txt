##
##	This file is part of qp42.
##
##	qp42 -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2012 by Janick Frasch, Hans Joachim Ferreau et al. 
##	All rights reserved.
##
##	qp42 is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qp42 is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qp42; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##

################################################################################
#
# Description:
#	CMake script for building qpDUNES static libraries, etc.
#
# Authors:
#	Milan Vukov, milan.vukov@esat.kuleuven.be
#
# Year:
#	2013.
#
# Usage:
#	- TBD
#
################################################################################


################################################################################
#
# Project settings
#
################################################################################

#
# Minimum required version of cmake 
#
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

#
# Project name and programming languages used
#
PROJECT( qpDUNES C CXX )

#
# Folder path for generated executables
#
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )

#
# Folder path for generated libraries
#
SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )

#
# CMake module(s) path
#
SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR} )

################################################################################
#
# Prerequisites
#
################################################################################

################################################################################
#
# User configuration
#
################################################################################

IF( NOT CMAKE_VERBOSE_MAKEFILE )
	SET( CMAKE_VERBOSE_MAKEFILE OFF )
ENDIF( NOT CMAKE_VERBOSE_MAKEFILE )

#
# Build type
#
IF( NOT CMAKE_BUILD_TYPE )
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
		FORCE
	)
ENDIF( NOT CMAKE_BUILD_TYPE )

OPTION( QPDUNES_WITH_LAPACK
	"Build qpOASES using original LAPACK routines"
	OFF
)

OPTION( QPDUNES_WITH_BLAS
	"Build qpOASES using original BLAS routines"
	OFF
)

OPTION( QPDUNES_WITH_OPENMP
	"Build qpDUNES with OpenMP support"
	OFF
)

################################################################################
#
# Compiler settings
#
################################################################################

SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall -Wextra -pedantic -Wshadow -g -finline-functions -DLINUX -D__NO_COPYRIGHT__ -D__DEBUG__ -D__FULL_LOGGING__ -std=c99" )

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wshadow -g -finline-functions -fPIC -DLINUX" )

IF ( QPDUNES_WITH_OPENMP )
	FIND_PACKAGE( OpenMP )
ENDIF()

IF( OPENMP_FOUND )
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ENDIF()

################################################################################
#
# Include directories
#
################################################################################

INCLUDE_DIRECTORIES(
	.
	./include
	./include/qp
	./interfaces
	./interfaces/mpc
	
	# qpOASES include paths
	#./externals/qpOASES-3.0beta
	#./externals/qpOASES-3.0beta/include
)

################################################################################
#
# Sources and headers
#
################################################################################

SET( qpDUNES_HEADERS
	${CMAKE_CURRENT_SOURCE_DIR}/include/qpDUNES.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/qp/dual_qp.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/qp/matrix_vector.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/qp/setup_qp.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/qp/stage_qp_solver_clipping.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/qp/stage_qp_solver_qpoases.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/qp/types.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/qp/utils.h
)

SET( qpDUNES_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/src/dual_qp.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/stage_qp_solver_clipping.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/stage_qp_solver_qpoases.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/matrix_vector.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/setup_qp.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/utils.c
)

SET( QPOASES_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/Bounds.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/Constraints.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/Flipper.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/Indexlist.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/Matrices.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/MessageHandling.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/Options.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/QProblem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/QProblemB.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/SQProblem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/SubjectTo.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/Utils.cpp
)

FILE( GLOB QPOASES_HEADERS 
	  ${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/include/*.hpp )

# find LAPACK
IF ( QPOASES_WITH_LAPACK )
	ENABLE_LANGUAGE( C )
	ENABLE_LANGUAGE( Fortran )

	FIND_PACKAGE( LAPACK )

	IF ( NOT LAPACK_FOUND )
		SET( QPOASES_SOURCES ${QPOASES_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/LAPACKReplacement.cpp )
	ENDIF ( NOT LAPACK_FOUND )

ELSE ( QPOASES_WITH_LAPACK )
	SET( QPOASES_SOURCES ${QPOASES_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/LAPACKReplacement.cpp )
ENDIF ( QPOASES_WITH_LAPACK )

# find BLAS
IF ( QPOASES_WITH_BLAS )
	ENABLE_LANGUAGE( C )
	ENABLE_LANGUAGE( Fortran )

	FIND_PACKAGE( BLAS )

	IF ( NOT BLAS_FOUND )
		SET( QPOASES_SOURCES ${QPOASES_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/BLASReplacement.cpp )
	ENDIF ( NOT BLAS_FOUND )

ELSE ( QPOASES_WITH_BLAS )
	SET( QPOASES_SOURCES ${QPOASES_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/externals/qpOASES-3.0beta/src/BLASReplacement.cpp )
ENDIF ( QPOASES_WITH_BLAS )


################################################################################
#
# Build applications, libraries and examples
#
################################################################################

#
# Build the qpDUNES stand-alone library
#
ADD_LIBRARY(
	qpDUNES STATIC
	
	# Sources and headers
	${qpDUNES_HEADERS}
	${qpDUNES_SOURCES}
)

#ADD_LIBRARY(
#	qpDUNES_qpOASES STATIC
#
#	# qpOASES sources
#	${QPOASES_SOURCES}
#)

IF ( QPOASES_WITH_LAPACK AND LAPACK_FOUND )
	TARGET_LINK_LIBRARIES(
		qpDUNES_qpOASES
		${LAPACK_LIBRARIES}
	)
ENDIF()

IF ( QPOASES_WITH_BLAS AND BLAS_FOUND )
	TARGET_LINK_LIBRARIES(
		qpDUNES_qpOASES
		${BLAS_LIBRARIES}
	)
ENDIF()

TARGET_LINK_LIBRARIES(
	qpDUNES
	#qpDUNES_qpOASES
)

#
# Build the mpcDUNES library
#
ADD_LIBRARY(
	mpcDUNES STATIC
	
	# Sources and headers
	${CMAKE_CURRENT_SOURCE_DIR}/interfaces/mpc/setup_mpc.c
	${CMAKE_CURRENT_SOURCE_DIR}/interfaces/mpc/setup_mpc.h
)
TARGET_LINK_LIBRARIES(
	mpcDUNES
	qpDUNES
)

#
# Build the examples
# NOTE: Assumption is that all examples are in C and that one example
# is in exactly one file.
#

FILE( GLOB qpDUNES_EXAMPLES ${CMAKE_CURRENT_SOURCE_DIR}/examples/*.c )
FOREACH( EXAMPLE ${qpDUNES_EXAMPLES} )
		 GET_FILENAME_COMPONENT( EXEC_NAME ${EXAMPLE} NAME_WE )
		 ADD_EXECUTABLE( example_${EXEC_NAME} ${EXAMPLE} )
		 TARGET_LINK_LIBRARIES(
			example_${EXEC_NAME}
			mpcDUNES
			qpDUNES
			)
		 SET_TARGET_PROPERTIES( example_${EXEC_NAME} 
		 	PROPERTIES
				OUTPUT_NAME "${EXEC_NAME}"		   
		 )
ENDFOREACH()

################################################################################
#
# Configuration of the bash script
#
################################################################################

SET( qpDUNES_ENV_INCLUDE_DIRS
	${PROJECT_SOURCE_DIR}
	${PROJECT_SOURCE_DIR}/include
	${PROJECT_SOURCE_DIR}/include/qp
	${PROJECT_SOURCE_DIR}/interfaces
	${PROJECT_SOURCE_DIR}/interfaces/mpc
	#${PROJECT_SOURCE_DIR}/externals/qpOASES-3.0beta
	#${PROJECT_SOURCE_DIR}/externals/qpOASES-3.0beta/include
)

SET( qpDUNES_ENV_LIBRARY_DIRS
	${LIBRARY_OUTPUT_PATH}
)

SET( qpDUNES_ENV_BUILD_DIR
	${PROJECT_BINARY_DIR}
)

SET( qpDUNES_ENV_STATIC_LIBRARIES
	qpDUNES mpcDUNES
)

SET( qpDUNES_ENV_HEADERS
	${qpDUNES_HEADERS}
	#${QPOASES_HEADERS}
)

SET( qpDUNES_ENV_SOURCES
	${qpDUNES_SOURCES}
	#${QPOASES_SOURCES}
)

CONFIGURE_FILE(
	"${PROJECT_SOURCE_DIR}/cmake/qpDUNES_env.sh.in"
	"${PROJECT_BINARY_DIR}/qpDUNES_env.sh"
	@ONLY
)
