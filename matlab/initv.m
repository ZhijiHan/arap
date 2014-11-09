function [ V2 ] = initv( V, W, cid, C )
  % Tao Du
  % Nov 8, 2014
  
  % Given vertices, weight, the index and value of all the constrained
  % vertices, returns initial vertices from naive Laplacian equation.
  
  % The length of cid should match the # of rows in C.
  % V2 is a (# vertices) x 3 matrix.

  % Naive Laplacian equation:
  % min \|LV2 - LV\|, s.t, V2(cid, :) = C. where L = -W.
  
  % Get the number of vertices.
  vnum = size(V, 1);
  
  % Compute fid, which is the indices of all the free variables.
  fid = setdiff(1 : vnum, cid);
  
  % Compute L.
  L = -W;
  
  % Rewrite LV2 = L(:, fid) * V2(fid, :) + L(:, cid) * C.
  A = L(:, fid);
  B = L(:, cid);
  
  % LV2 - LV = A * V2(fid, :) - (LV - B * C)
  % Solve V2(fid, :) by \.
  V2(fid, :) = A \ (L * V - B * C);
  
  % Write back C into V2 as hard constraints.
  V2(cid, :) = C;
end

