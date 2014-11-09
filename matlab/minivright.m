function [ rhs ] = minivright( V, R, cid, C, N, W )
  % Tao Du
  % Nov 8, 2014
  
  % Given vertices, rotations, constrained vertices conditions,
  % neighborhoods and weights, try to compute the right hand side for 
  % minimizing new vertices.
  
  % Get the number of vertices.
  vnum = size(V, 1);
  
  % Compute fid, which is the indices of all the free variables.
  fid = setdiff(1 : vnum, cid);
  
  % The Laplacian equation.
  % For each i:
  % \sum_{j\in N(i)} w_{ij}(p_i' - p_j') = \sum_{j\in N(i)}
  % w_{ij} / 2 * (R_i + R_j).
  
  % Initialize the right hand side with fixed vertices
  rhs = W(fid, cid) * C;
  
  % Loop over all the free variables.
  % Set base to record the row we want to update rhs.
  base = 1;
  for i = fid
    % Get all the vertices incident on i.
    incidence = find(N(i, :));
    
    % Loop over all i's neighbors.
    for j = incidence
      % Get weight.
      w = W(i, j);
      
      % Get rotation matrix.
      basei = (i - 1) * 3;
      ri = R(basei + 1 : basei + 3, :);
      basej = (j - 1) * 3;
      rj = R(basej + 1 : basej + 3, :);
      
      % Update the base-th row in rhs.
      rhs(base, :) = rhs(base, :) + ...
                     w / 2 * (V(i, :) - V(j, :)) * (ri + rj)';
    end
    base = base + 1;
  end
end

