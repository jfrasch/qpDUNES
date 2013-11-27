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
%  A simple double integrator model for a badminton robot
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

INFTY = 1.0e12;
QUADPROG = true;
dt = 0.01;              % discretization sampling time 10ms, needed for constraints
	

% Problem data


Q = [[ 1.e-4, 0.0   ]
	 [ 0.0  , 1.e-4 ]];
R = [[ 1.0e0 ]];
P = Q;
	
A = [[ 1.0, 1.0*dt ]
	 [ 0.0, 1.0    ]];	
B = [[ 0.0    ]
     [ 1.0*dt ]];
c = [ 0.0;
      0.0 ];

xLow = [ -1.9, -3.0 ]';
xUpp = [  1.9,  3.0 ]';
uLow = [ -30.0 ]';
uUpp = [  30.0 ]';
xRef = [  0.0,  0.0 ]';
uRef = [  0.0 ]';


% build up data
Hi = blkdiag(Q, R);
H = repmat( Hi, 1, nI );
Ci = [A B];
C = repmat( Ci, 1, nI );
cFull = repmat( c, nI, 1 );

ziLow = [ xLow; uLow ];
ziUpp = [ xUpp; uUpp ];
ziRef = [ xRef; uRef ];
zLow = [ repmat( ziLow, nI, 1 ); xLow ];
zUpp = [ repmat( ziUpp, nI, 1 ); xUpp ];
zRef = [ repmat( ziRef, nI, 1 ); xRef ];


% arrival constraints
idxArrivalStart = 44;       			% 43 and shorter is infeasible
idxArrivalEnd = idxArrivalStart + 1;

if (idxArrivalEnd < nI)
    k = idxArrivalStart;
    zLow(k*nZ+1) = zRef(k*nZ+1);
    zUpp(k*nZ+1) = zRef(k*nZ+1);
    k = idxArrivalEnd;
    zLow(k*nZ+1) = zRef(k*nZ+1);
    zUpp(k*nZ+1) = zRef(k*nZ+1);
end


% initial value and QP solution

x0 = [ -1, 0 ];


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

	
disp( ['Solving double integrator [nI = ', num2str(nI), ', nX = ', num2str(nX), ', nU = ', num2str(nU), ']'] );


% -- CLEANUP --
mpcDUNES( 'cleanup' );		% not really needed, just to be safe


% -- INIT z-STYLE --
mpcDUNES( 'init', nI, ...
          H, P, [], C, cFull, ...
          zLow, zUpp, zRef, ...
          qpOptions );
% -- UPDATE --
% mpcDUNES( 'update', ...
%           H, P, [], C, c, ...
%           [], zLow, zUpp, zRef );
% -- SOLVE --
[uOpt, xOpt, stat, objFctnVal] = mpcDUNES( 'solve', x0 );
% -- CLEANUP --
mpcDUNES( 'cleanup' );


%% plotting

t = dt * (1:nI+1);
linewidth = 2;

figure; 
subplot(3,1,1); hold on;
plot(t,xOpt(1:2:end),'b', 'linewidth', linewidth); 
plot([t(1) t(end)],[xLow(1) xLow(1)], 'm-', 'linewidth', linewidth);     % lb
plot([t(1) t(end)],[xUpp(1) xUpp(1)], 'm-', 'linewidth', linewidth);     % ub
subplot(3,1,2); hold on;
plot(t,xOpt(2:2:end),'b', 'linewidth', linewidth); 
plot([t(1) t(end)],[xLow(2) xLow(2)], 'm-', 'linewidth', linewidth);     % lb
plot([t(1) t(end)],[xUpp(2) xUpp(2)], 'm-', 'linewidth', linewidth);     % ub
subplot(3,1,3); hold on;
plot(t(1:end-1),uOpt(1:end),'b', 'linewidth', linewidth);
plot([t(1) t(end)],[uLow(1) uLow(1)], 'm-', 'linewidth', linewidth);     % lb
plot([t(1) t(end)],[uUpp(1) uUpp(1)], 'm-', 'linewidth', linewidth);     % ub

% Arrival constraint
subplot(3,1,1); 

pos = [t(idxArrivalStart) ...                   % x
       zUpp(nZ*idxArrivalStart+1) ...           % y
       t(idxArrivalEnd)-t(idxArrivalStart) ...  % w
       zUpp(nZ*(idxArrivalStart-1)+1)-zUpp(nZ*idxArrivalStart+1) ...    % h
       ];
rectangle('pos',pos ,'facecolor',[.8 .8 .8],'edgecolor','none');
pos = [t(idxArrivalStart) ...                   % x
       zLow(nZ*(idxArrivalStart-1)+1) ...       % y
       t(idxArrivalEnd)-t(idxArrivalStart) ...  % w
       zLow(nZ*idxArrivalStart+1)-zLow(nZ*(idxArrivalStart-1)+1) ...    % h
       ];
rectangle('pos',pos ,'facecolor',[.8 .8 .8],'edgecolor','none');
 
hold off;


%% end of file



 
 