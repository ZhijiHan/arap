function [ N ] = neighbor( F )
  % Tao Du
  % Nov 8, 2014
  
  % Given a face matrix, compute an adjacent matrix for neighborhood. The
  % return matrix N is a sparse logical matrix, with N(i, j) = 1 indicates
  % there is an edge between vertex i and j. N is guaranteed to be
  % symmetric.
  
  % Get the largest index in F.
  vnum = max(max(F));
  
  % Get indices from three columns in F.
  X = F(:, 1);
  Y = F(:, 2);
  Z = F(:, 3);
  
  % Generate neighbor matrix.
  N = sparse([X; Y; X; Z; Y; Z], [Y; X; Z; X; Z; Y], 1, vnum, vnum);
  
  % transform N into a logic matrix.
  N = N > 0;
end

