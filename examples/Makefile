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



##
##	Filename:  examples/Makefile
##	Author:    Janick Frasch, Hans Joachim Ferreau
##	Version:   1.0beta
##	Date:      2012
##

SRCDIR = ../src
INTERFACEDIR = ../interfaces

# select your operating system here!
include ../make_linux.mk
#include ../make_windows.mk


##
##	flags
##

IFLAGS      =  -I. \
               -I../include \
               -I../interfaces

QPDUNES_EXES = \
	example1${EXE}



##
##	targets
##

all: ${QPDUNES_EXES}


example1${EXE}: example1.${OBJEXT} ../src/libqpdunes.a
	${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPDUNES_LIB} ${LIBS}




clean:
	${RM} -f *.${OBJEXT} ${QP42_EXES}

clobber: clean


%.${OBJEXT}: %.c
	@echo "Creating" $@
	${CC} ${DEF_TARGET} ${IFLAGS} ${CCFLAGS} -c $<


%.${OBJEXT}: %.cpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} ${IFLAGS} ${CPPFLAGS} -c $<


##
##	end of file
##
