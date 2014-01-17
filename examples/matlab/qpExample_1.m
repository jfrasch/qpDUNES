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

%% Abstract
%
%  A simple QP example
%
%  execute
%     cd ../../interfaces/matlab/; make; cd ../../examples/matlab
%  to compile mex
%


%% Problem definition
% Problem dimensions

nI = 80;            % number of control intervals
nX = 2;             % number of states
nU = 1;             % number of controls
nZ = nX+nU;
nD = zeros(nI+1,1); % number of affine constraints



%  Problem settings

dt = 0.01;              % discretization sampling time 10ms, needed for constraints
	

% Problem data

Hi = [[ 1.e-4, 0.0  , 0.0   ]
      [ 0.0  , 1.e-4, 0.0   ]
      [ 0.0  , 0.0  , 1.0e0 ]];
gi = [ 0.0;
       0.0;
       0.0 ];
P = [[ 1.e-4, 0.0 ]
     [ 0.0  , 1.e-4    ]];
	
Ci = [[ 1.0, 1.0*dt, 0.0    ]
      [ 0.0, 1.0   , 1.0*dt ]];	
ci = [ 0.0;
       0.0 ];

ziLow = [ -1.9, -3.0, -30.0 ]';
ziUpp = [  1.9,  3.0,  30.0 ]';


% build up data
H = repmat( Hi, 1, nI );
g = [ repmat( gi, nI, 1 ); gi(1:nX) ];
C = repmat( Ci, 1, nI );
c = repmat( ci, nI, 1 );

zLow = [ repmat( ziLow, nI, 1 ); ziLow(1:nX) ];
zUpp = [ repmat( ziUpp, nI, 1 ); ziUpp(1:nX) ];


% initial and terminal values
x0 = [ -1, 0 ];
xf = [ 0, 0 ];


%% SOLVE
% Options

% add qpDUNES path
qpDUNES_PATH = '../../interfaces/matlab';
addpath(qpDUNES_PATH);

qpOptions = qpDUNES_options( 'default', ...
                             'maxIter', 100, ...
                             'printLevel', 2, ...
                             'logLevel', 0, ...     % log all data
                             'lsType', 4, ...       % Accelerated gradient biscection LS
                             'stationarityTolerance', 1.e-6, ...
                             'regType', 2 ...       % regularize only singular directions; 1 is normalized Levenberg Marquardt
                             ...
                             );


% a) initialize data
qpDUNES( 'init', nI, ...
          H, P, g, ...
          C, c, ...
          zLow, zUpp, ...
          [], [], [], ...        % affine constraints not yet supported
          qpOptions );
% b) now add an initial value constraint (update on stage index 0)
qpDUNES( 'stageUpdate', 0, [], [], [], [], [x0,ziLow(nX+1:nZ)], [x0,ziUpp(nX+1:nZ)], [], [], [] );
% c) solve again
[zOpt, stat, lambda, mu, objFctnVal] = qpDUNES( 'solve' );

% d) add a terminal constraint (update on stage index nI)
qpDUNES( 'stageUpdate', nI, [], [], [], [], xf, xf, [], [], [] );
% e) solve yet once more
[zOpt, stat, lambda, mu, objFctnVal] = qpDUNES( 'solve' );

% % f) optional: do a shift and a full data update
% qpDUNES( 'shift');
% qpDUNES( 'update', H, P, g, C, c, zLow, zUpp, [], [], [] );
% qpDUNES( 'stageUpdate', 0, [], [], [], [], [x0,ziLow(nX+1:nZ)], [x0,ziUpp(nX+1:nZ)], [], [], [] );
% [zOpt, stat, lambda, mu, objFctnVal] = qpDUNES( 'solve' );


%% plotting

t = dt * (1:nI+1);
linewidth = 1;

figure; 
subplot(3,1,1); plot(t,zOpt(1:3:end),'b', 'linewidth', linewidth); 
subplot(3,1,2); plot(t,zOpt(2:3:end),'b', 'linewidth', linewidth); 
subplot(3,1,3); stairs(t(1:end-1),zOpt(3:3:end),'r', 'linewidth', linewidth);

hold off;


%% end of file



 
 