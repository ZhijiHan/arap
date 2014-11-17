% Tao Du
% Nov 8, 2014

% This script implements the arap algorithm.

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
R = minir(V, V2, N, W);

% Compute arap energy.
arap = comparap(V, V2, R, N, W);
disp(arap);

% Build a temporary folder for displaying.
showmodel(V, F, C);

% Set the max number of iteration.
maxiter = 50;
iter = 0;

% Prefactorize the left matrix for minimizing vertices.
L = minivleft(W, cid);
% LU factorization.
[l, u] = lu(L);

while iter < maxiter
  % Solve rotations.
  R = minir(V, V2, N, W);

  % Compute the right hand side for minimizing vertices.
  rhs = minivright(V, R, cid, C, N, W);

  % Write free solutions back to V2.
  % lux = rhs.
  % ly = rhs.
  % ux = y.
  y = l \ rhs;
  V2(fid, :) = u \ y;

  % Compute arap energy.
  arap = comparap(V, V2, R, N, W);
  disp(arap);

  % Display the model.
  showmodel(V2, F, C);

  % Increment iter.
  iter = iter + 1;
end

clear all;