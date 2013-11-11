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
##	Filename:  Makefile
##	Author:    Janick Frasch, Hans Joachim Ferreau
##	Version:   1.0beta
##	Date:      2012
##


##
##	targets
##


all:
	@  cd src               			&& ${MAKE} && cd .. \
	&& cd interfaces/mpc    			&& ${MAKE} && cd ../.. \
	&& cd examples          			&& ${MAKE} && cd .. 

clean:
	@  cd src               			&& ${MAKE} clean && cd .. \
	&& cd interfaces/mpc    			&& ${MAKE} clean && cd ../.. \
	&& cd examples          			&& ${MAKE} clean && cd .. 

clobber: clean


.PHONY : all clean clobber doc testing


##
##   end of file
##
