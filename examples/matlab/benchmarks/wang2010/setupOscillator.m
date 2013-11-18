%	This file is part of qpDUNES.
%
%	qpDUNES -- An Implementation of the Dual Nonsmooth Newton Strategy.
%	Copyright (C) 2012 by Janick Frasch, Hans Joachim Ferreau et al. 
%	All rights reserved.
%
%	qpDUNES is free software; you can redistribute it and/or
%	modify it under the terms of the GNU Lesser General Public
%	License as published by the Free Software Foundation; either
%	version 2.1 of the License, or (at your option) any later version.
%
%	qpDUNES is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%	See the GNU Lesser General Public License for more details.
%
%	You should have received a copy of the GNU Lesser General Public
%	License along with qpDUNES; if not, write to the Free Software
%	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%


%% Abstract
%
%  As in Wang2010, Ferreau2012b
%


k = 1;          % spring constant 
lam = 0;        % damping constant 
a = -2*k;
b = -2*lam;
c = k;
d = lam;

n = 12; % state dimension
m = 3; % input dimension

Acts = [zeros(n/2) eye(n/2);
    [a,c,0,0,0,0,b,d,0,0,0,0;
     c,a,c,0,0,0,d,b,d,0,0,0;
     0,c,a,c,0,0,0,d,b,d,0,0;
     0,0,c,a,c,0,0,0,d,b,d,0;
     0,0,0,c,a,c,0,0,0,d,b,d;
     0,0,0,0,c,a,0,0,0,0,d,b]];

Bcts = [zeros(n/2,m);
    [1, 0, 0;
    -1, 0, 0;
     0, 1, 0;
     0, 0, 1;
     0,-1, 0;
     0, 0,-1]];

% convert to discrete-time system
ts = 0.5;       % sampling time
A = expm(ts*Acts);
B = (Acts\(expm(ts*Acts)-eye(n)))*Bcts;

% objective matrices
Q = eye(n);      
R = eye(m);        

% state and control limits
Xmax = 4; Umax = 0.5;
xmin = -Xmax*ones(n,1);
xmax = Xmax*ones(n,1);
umin = -Umax*ones(m,1);
umax = Umax*ones(m,1);
