function [ R ] = minir( V, V2, N, W )
  % Tao Du
  % Nov 8, 2014
  
  % Given vertices, new vertices, neighborhood and weight, minimize
  % rotation matrices by SVD projection.
  
  % The return value R is a (#vertices x 3) x 3 matrix.
  
  % Get the number of vertices.
  vnum = size(V, 1);
  
  % Get all the neighbors.
  [I, J, ~] = find(N);
  
  % Preallocate space for temperary matrix used to solve rotations.
  E = zeros(3 * vnum, 3);
  
  % Get the number of edges.
  n = length(I);
  
  % Loop over all the neighborhoods.
  for e = 1 : n
    % Get the index of the two incident vertices.
    i = I(e);
    j = J(e);
    
    % Get the old and new vertices.
    vni = V2(i, :)';
    vnj = V2(j, :)';
    voi = V(i, :)';
    voj = V(j, :)';
    
    % Get weight.
    w = W(i, j);
    
    % Compute the base index of the rotation matrix in E.
    base = (i - 1) * 3;
    
    % Add the contribution of this edge into E.
    E(base + 1 : base + 3, :) = E(base + 1 : base + 3, :) + ...
                                w * (voi - voj) * (vni - vnj)';
  end
  
  % Preallocate space for rotation matrix.
  R = zeros(3 * vnum, 3);
  
  % Use SVD to solve rotation matrix.
  for i = 1 : vnum
    % Compute the base index for vertex i.
    base = (i - 1) * 3;
    [u, ~, v] = svd(E(base + 1 : base + 3, :));
    r = v * u';
    
    % If det(r) is negative, we flip the last column in u.
    if det(r) < 0
      u(:, end) = -1 * u(:, end);
      r = v * u';
    end
    
    % Write r back to R.
    R(base + 1 : base + 3, :) = r;
  end
end

