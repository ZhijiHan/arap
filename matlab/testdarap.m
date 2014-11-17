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
R = minir(V, V2, N, W);

% Initialize V3.
V3 = V2(fid, :);

% Call comparapb.
[arap, darap] = comparapb(V3, V, cid, C, N, W);
disp('arap:');
disp(arap);

% Reshape darap because it is now a vector.
darap = reshape(darap, length(darap) / 3, 3);

% Now, test with darap with forward-backward check.
epsilon = 1e-8;
tolerance = 0.03;

% Collect all the relative errors.
rels = zeros(length(fid), 3);
for i = 1 : length(fid)
  % Loop over all the free variables.
  for j = 1 : 3
    % Cache V2(fid, j)
    cache = V2(fid(i), j);
    
    % Increment it.
    V2(fid(i), j) = cache + epsilon;
    arap2 = comparap(V, V2, R, N, W);
    
    % Decrement it.
    V2(fid(i), j) = cache - epsilon;
    arap3 = comparap(V, V2, R, N, W);
    
    % Reset.
    V2(fid(i), j) = cache;
    
    % Compute numerical derivatives.
    numer = (arap2 - arap3) / (2 * epsilon);
    
    % Get analytical derivatives.
    analy = darap(i, j);
    
    disp('Numerical derivatives:');
    disp(numer);
    disp('Analytic derivatives:');
    disp(analy);
    
    % Compute relative error.
    if analy ~= 0
      disp('Relative error:');
      rel = abs((numer - analy) / analy);
      disp(rel);
      rels(i, j) = rel;
      
      % If relative error is greater than 1%, complain
      if rel > tolerance
        disp('Test failed!');
      end
    end
  end
end

if max(max(rels)) < tolerance
  disp('All tests passed!');
end
