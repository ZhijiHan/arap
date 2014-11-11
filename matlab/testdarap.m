% Tao Du
% Nov 11, 2014

% This script tries to test whether our implementation for comparapb is
% correct.

% Clear variables, reset display format to long.
clc; clear all;
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

% Initialize V3.
V3 = V2(fid, :);

% Call comparapb.
[arap, darap] = comparapb(V3, V, cid, C, N, W);
disp('arap:');
disp(arap);

clear all;