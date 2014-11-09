function [ L ] = minivleft( W, cid )
  % Tao Du
  % Nov 8, 2014
  
  % This function computes the left matrix for minimizing vectors.
  
  % Get the # vertices.
  vnum = size(W, 1);
  
  % Get index of the free variables.
  fid = setdiff(1 : vnum, cid);
  
  % Compute L.
  L = -W(fid, fid);
end

