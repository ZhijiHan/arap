function [ W ] = compweight( V, F )
  % Tao Du
  % Nov 8, 2014
  
  % Computes the cotangent weight for given V and F.
  % The output W is a sparse matrix. W(i, j) is the cotangent weight for
  % edge (i, j) in the triangle mesh. W is guaranteed to be symmetric.
  % edge(i, i) is not defined. We use edge(i, i) to store the negative sum
  % of all the edge(i, j), i.e. edge(i, i) = -\sum_{j \neq i} edge(i, j)
  
  % Get the vertex number.
  vnum = size(V, 1);
  
  % Get the face number;
  fnum = size(F, 1);
  
  % Preallocate index i, j and w.
  i = zeros(12 * fnum, 1);
  j = zeros(12 * fnum, 1);
  w = zeros(12 * fnum, 1);
  
  % Loop over all the triangles.
  for f = 1 : fnum
    % Extract the index of three vertices.
    A = F(f, 1);
    B = F(f, 2);
    C = F(f, 3);
    
    % Extract the triangles and half cotangent values. 
    T = V(F(f, :), :);
    half = tricot(T) / 2;
    
    % Insert new weights for edge BC. We should add weights for edge BC,
    % edge CB, vertex B and vertex C.
    base = (f - 1) * 12;
    i(base + 1 : base + 4) = [B; C; B; C];
    j(base + 1 : base + 4) = [C; B; B; C];
    w(base + 1 : base + 4) = [half(1); half(1); -half(1); -half(1)];
    
    % Insert new weights for edge AC.
    i(base + 5 : base + 8) = [A; C; A; C];
    j(base + 5 : base + 8) = [C; A; A; C];
    w(base + 5 : base + 8) = [half(2); half(2); -half(2); -half(2)];
    
    % Insert new weights for edge AB.
    i(base + 9 : base + 12) = [A; B; A; B];
    j(base + 9 : base + 12) = [B; A; A; B];
    w(base + 9 : base + 12) = [half(3); half(3); -half(3); -half(3)];
  end

  % Build sparse weight matrix.
  W = sparse(i, j, w, vnum, vnum);
end

