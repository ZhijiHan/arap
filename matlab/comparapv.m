function [ arap ] = comparapv( V, V2, N, W )
  % Tao Du
  % Nov 8, 2014
  
  % Given the old and new vertices, neighborhood and weight, computes the 
  % arap energy.
  % The rotation matrix is given by minir.
  
  % Compute R.
  R = minir(V, V2, N, W);
  
  % Get the nonzero edges.
  [I, J, ~] = find(N);
  
  % Get the number of edges.
  enum = length(I);
  
  % Initialize arap energy.
  arap = 0.0;
  
  % Loop over all the edges.
  for e = 1 : enum
    % Get the indices of the first and second edges.
    i = I(e);
    j = J(e);
    
    % Get the weight of the edge.
    w = W(i, j);
    
    % Get the old and new vertices.
    voi = V(i, :);
    voj = V(j, :);
    vni = V2(i, :);
    vnj = V2(j, :);
    
    % Get the rotation matrix.
    base = (i - 1) * 3;
    r = R(base + 1 : base + 3, :);
    
    % Add this edge into arap energy.
    arap = arap + compedge(w, vni, vnj, r, voi, voj);
  end
end

