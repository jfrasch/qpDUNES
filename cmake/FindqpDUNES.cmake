##
##	This file is part of qpDUNES.
##
##	qpDUNES -- A DUal NEwton Strategy for convex quadratic programming.
##	Copyright (C) 2012-2014 by Janick Frasch, Hans Joachim Ferreau et al. 
##	All rights reserved.
##
##	qpDUNES is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qpDUNES is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qpDUNES; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##


################################################################################
#
# Description:
#	qpDUNES package configuration file
#
#
# Authors:
#	Milan Vukov, milan.vukov@esat.kuleuven.be
#
# Year:
#	2013 - 2014.
#
# NOTE:
#	- This script is for Linux/Unix use only.
#
#	- PREREQUISITE: sourced qpDUNES_env.sh in your ~/.bashrc file. This script
#		will try to find qpDUNES folders, libraries etc., but looking for them
#		in environmental variables.
#
# Usage:
#	- Linux/Unix: TODO
#
################################################################################

################################################################################
#
# Search for package components
#
################################################################################

MESSAGE( STATUS "********************************************************************************" )
MESSAGE( STATUS "Looking for qpDUNES package: \n" )

#
# Include folders
#
MESSAGE( STATUS "Looking for qpDUNES include directories" )
SET( qpDUNES_INCLUDE_DIRS $ENV{qpDUNES_ENV_INCLUDE_DIRS} )
IF( qpDUNES_INCLUDE_DIRS )
	MESSAGE( STATUS "Found qpDUNES include directories: ${qpDUNES_INCLUDE_DIRS} \n" )
	SET( qpDUNES_INCLUDE_DIRS_FOUND TRUE )
ELSE( qpDUNES_INCLUDE_DIRS )
	MESSAGE( STATUS "Could not find qpDUNES include directories \n" )
ENDIF( qpDUNES_INCLUDE_DIRS )

#
# Library folders
#
MESSAGE( STATUS "Looking for qpDUNES library directories" )
SET( qpDUNES_LIBRARY_DIRS $ENV{qpDUNES_ENV_LIBRARY_DIRS} )
IF( qpDUNES_LIBRARY_DIRS )
	MESSAGE( STATUS "Found qpDUNES library directories: ${qpDUNES_LIBRARY_DIRS} \n" )
	SET( qpDUNES_LIBRARY_DIRS_FOUND TRUE )
ELSE( qpDUNES_LIBRARY_DIRS )
	MESSAGE( STATUS "Could not find qpDUNES library directories \n" )
ENDIF( qpDUNES_LIBRARY_DIRS )

#
# Libraries
#
FIND_LIBRARY( qpDUNES_STATIC_LIBRARIES
	NAMES qpdunes
	PATHS ${qpDUNES_LIBRARY_DIRS}
	NO_DEFAULT_PATH
)
IF( qpDUNES_STATIC_LIBRARIES )
	MESSAGE( STATUS "Found qpDUNES static library: ${qpDUNES_STATIC_LIBRARIES} \n" )
ELSE( qpDUNES_STATIC_LIBRARIES )
	MESSAGE( STATUS "Could not find qpDUNES static library.\n" )
	SET( qpDUNES_STATIC_LIBS_FOUND FALSE )
ENDIF( qpDUNES_STATIC_LIBRARIES )

#
# And finally set found flag...
#
IF( qpDUNES_INCLUDE_DIRS_FOUND AND qpDUNES_LIBRARY_DIRS_FOUND AND qpDUNES_STATIC_LIBS_FOUND )
	SET( qpDUNES_FOUND TRUE )
ENDIF()

MESSAGE( STATUS "********************************************************************************" )