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
##	Filename:  src/make_linux
##	Author:    Janick Frasch, Hans Joachim Ferreau
##	Version:   1.0beta
##	Date:      2012
##


##
##	definitions for compiling with gcc under linux
##

CC = gcc
CPP = g++
AR  = ar
RM  = rm

OBJEXT = o
LIBEXT = a
EXE =
DEF_TARGET = -o $@

OMPFLAGS = 


CCFLAGS = -Wall -pedantic -Wshadow -O3 -finline-functions -DLINUX -U__MEASURE_TIMINGS__ -U__ANALYZE_FACTORIZATION__ -std=c99		##C99 temporary to avoid warnings

QPDUNES_LIB         =  -L${SRCDIR} -lqpdunes

MPCDUNES_LIB        =  -L${INTERFACEDIR}/mpc -lmpcDUNES


LIBS         =  -lm


##
##	end of file
##
