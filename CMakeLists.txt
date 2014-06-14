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
#	2013 - 2014.
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
PROJECT( qpDUNES C )

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
# User configuration
#
################################################################################

#
# Build type
#
IF( NOT CMAKE_BUILD_TYPE )
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
		FORCE
	)
ENDIF( NOT CMAKE_BUILD_TYPE )

OPTION( QPDUNES_WITH_OPENMP
	"Build qpDUNES with OpenMP support"
	OFF
)

################################################################################
#
# Compiler settings
#
################################################################################

SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall -Wextra -pedantic -Wshadow -g -finline-functions -DLINUX -D__NO_COPYRIGHT__ -D__DEBUG__ -D__FULL_LOGGING__ -D__SIMPLE_BOUNDS_ONLY__ -std=c99" )
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")

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
	# mpcDUNES
	${CMAKE_CURRENT_SOURCE_DIR}/interfaces/mpc/setup_mpc.h
)

SET( qpDUNES_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/src/dual_qp.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/stage_qp_solver_clipping.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/stage_qp_solver_qpoases.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/matrix_vector.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/setup_qp.c
	${CMAKE_CURRENT_SOURCE_DIR}/src/utils.c
	# mpcDUNES
	${CMAKE_CURRENT_SOURCE_DIR}/interfaces/mpc/setup_mpc.c
)

################################################################################
#
# Build applications, libraries and examples
#
################################################################################

#
# Build the qpDUNES stand-alone library
#
ADD_LIBRARY(
	qpdunes STATIC
	
	# Sources and headers
	${qpDUNES_HEADERS}
	${qpDUNES_SOURCES}
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
		qpdunes
	)
	SET_TARGET_PROPERTIES( example_${EXEC_NAME} 
		PROPERTIES
		OUTPUT_NAME "${EXEC_NAME}"		   
	)
ENDFOREACH()

################################################################################
#
# Installation rules
#
################################################################################

INSTALL(
    DIRECTORY
	    include/
    DESTINATION
        include/qpdunes
)

INSTALL(
    FILES
	    interfaces/mpc/setup_mpc.h
    DESTINATION
        include/qpdunes/interfaces/mpc
)

INSTALL(
	TARGETS
	    qpdunes
    DESTINATION
        lib
)

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
)

SET( qpDUNES_ENV_LIBRARY_DIRS
	${LIBRARY_OUTPUT_PATH}
)

SET( qpDUNES_ENV_BUILD_DIR
	${PROJECT_BINARY_DIR}
)

SET( qpDUNES_ENV_STATIC_LIBRARIES
	qpdunes
)

SET( qpDUNES_ENV_HEADERS
	${qpDUNES_HEADERS}
)

SET( qpDUNES_ENV_SOURCES
	${qpDUNES_SOURCES}
)

CONFIGURE_FILE(
	"${PROJECT_SOURCE_DIR}/cmake/qpDUNES_env.sh.in"
	"${PROJECT_BINARY_DIR}/qpDUNES_env.sh"
	@ONLY
)
