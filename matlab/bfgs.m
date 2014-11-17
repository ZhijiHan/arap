% Tao Du
% Nov 12, 2014

% This script implements the bfgs algorithm.

% Clear variables, reset display format to long.
clear all;
format long;

% Read square_21_spikes.
[V, F] = readoff('../model/square_21_spikes.off');

% Get vertex number.
vnum = size(V, 1);

% Read cid and constrained vertices.
id = readdmat('../model/square_21_spikes-selection.dmat');
index = 1 : vnum;
cid = index(id > 0);
C = readdmat('../model/square_21_spikes.dmat');

% Compute the indices of free variables.
fid = setdiff(1 : vnum, cid);

% Compute weight.
W = compweight(V, F);

% Compute neighbors.
N = neighbor(F);

% Initialize V and R.
V2 = initv(V, W, cid, C);
R = minir(V, V2, N, W);
V3 = V2(fid, :);

% Compute arap energy.
arap = comparap(V, V2, R, N, W);
disp(arap);

% Set the max number of iteration.
maxiter = 50;

% Calling minFunc.
options = [];
options.Display = 'full';
options.MaxFunEvals = maxiter;
options.Method = 'lbfgs';
options.DerivativeCheck = 'off';

% After maxiter iterations, get the optimal solution.
[V3, arap] = minFunc(@comparapb, V3(:), options, V, cid, C, N, W);
V2(fid, :) = reshape(V3, length(V3) / 3, 3);
disp(arap);