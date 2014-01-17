%%
%	This file is part of qpDUNES.
%
%	qpDUNES -- A DUal NEwton Strategy for convex quadratic programming.
%	Copyright (C) 2012 by Janick Frasch, Hans Joachim Ferreau et al. 
%	All rights reserved.
%
%	qpDUNES is free software; you can redistribute it and/or
%	modify it under the terms of the GNU Lesser General Public
%	License as published by the Free Software Foundation; either
%	version 2.1 of the License, or (at your option) any later version.
% 
% 	qpDUNES is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% 	See the GNU Lesser General Public License for more details.
% 
% 	You should have received a copy of the GNU Lesser General Public
% 	License along with qpDUNES; if not, write to the Free Software
% 	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
 



%
%	Filename:  interfaces/matlab/make.m
%	Author:    Janick Frasch, Hans Joachim Ferreau
%	Version:   1.0beta
%	Date:      2012
%


%%
% consistency check
if ( exist( [pwd, '/make.m'],'file' ) == 0 )
	disp( 'ERROR: Run this make script directly within the directory' );
	disp( '       <qpDUNES-inst-dir>/interfaces/matlab, please.' );
	return;
end


QPDUNESPATH = '../../';

IFLAGS  = [ '-I.',' ', ...
            '-I',QPDUNESPATH,'include',' ', ...
            '-I',QPDUNESPATH,'externals/qpOASES-3.0beta/include',' ', ...
            '-I',QPDUNESPATH,'interfaces/mpc',' ' ];


if ( ispc == 0 )
	CPPFLAGS  = [ IFLAGS, '-largeArrayDims -D__cpluplus -D__MATLAB__ -cxx -O -DLINUX CFLAGS=''$CFLAGS -fPIC''', ' ' ]; 
else
	CPPFLAGS  = [ IFLAGS, '-largeArrayDims -D__cpluplus -D__MATLAB__ -O -DWIN32', ' ' ];
end

QPDUNES_OBJECTS =	[	QPDUNESPATH, 'src/stage_qp_solver_clipping.c ',...
						QPDUNESPATH, 'src/stage_qp_solver_qpoases.cpp ',...
						QPDUNESPATH, 'src/utils.c ',...
						QPDUNESPATH, 'src/dual_qp.c ',...
						QPDUNESPATH, 'src/matrix_vector.c ',...
						QPDUNESPATH, 'src/setup_qp.c ',...
						QPDUNESPATH, 'interfaces/mpc/setup_mpc.c ',...
						];

DEBUGFLAGS = ' ';
% DEBUGFLAGS = ' -g -D__DEBUG__ -v CXXDEBUGFLAGS=''$CXXDEBUGFLAGS -D__DEBUG__ -Wall -pedantic -Wshadow'' ';

MEXOBJS = { 'qpDUNES', ...
			'mpcDUNES' 
			}

;
for i = 1:length(MEXOBJS)
    cmd = [ 'mex -v CC="gcc" CXX="gcc" LD="gcc" COPTIMFLAGS="$COPTIMFLAGS -O3" -output ', MEXOBJS{i}, ' ', DEBUGFLAGS, CPPFLAGS, [MEXOBJS{i},'.cpp ',QPDUNES_OBJECTS] ]
	
	eval( cmd );
	disp( [ MEXOBJS{i}, '.', eval('mexext'), ' successfully created!'] );
end

path( path,pwd );


clear QPDUNESPATH IFLAGS CPPFLAGS QPDUNES_OBJECTS DEBUGFLAGS MEXOBJS cmd



%%
%%	end of file
%%
