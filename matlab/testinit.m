% Tao Du
% Nov 8, 2014

% Test initv and minir implementations.

% Read square_21_spikes.
clc; clear all;
[V, F] = readoff('../model/square_21_spikes.off');

% Get vertex number.
vnum = size(V, 1);

% Read cid and constrained vertices.
id = readdmat('../model/square_21_spikes-selection.dmat');
index = 1 : vnum;
cid = index(id > 0);
C = readdmat('../model/square_21_spikes.dmat');

% Compute weight.
W = compweight(V, F);

% Compute neighbors.
N = neighbor(F);

% Initialize V and R.
V2 = initv(V, W, cid, C);
R = minir(V, V2, N, W);

% Compute arap energy.
arap = comparap(V, V2, R, N, W);

% Display the arap energy.
disp('arap = ');
disp(arap);
clear all;