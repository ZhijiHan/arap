function [ arap, darap ] = comparapb( V3, V, cid, C, N, W )
  % Tao Du
  % Nov 9, 2014
  
  % This function is used in bfgs optimization. Unfortunately we have to
  % follow the convention that unknown variables are the first parameters.
  % The return values arap is the value of the arap energy, and darap is
  % the gradient of arap w.r.t all the vertices.
  
  % Note that now the definition of arap changes because we treat R as a
  % function of vertices.
  % arap = \sum w_{ij} \|p_i' - p_j' - r_i(p)(p_i - p_j)\|^2
  
  % V3 contains only the unknowns vertices. You can think of V3 as V3 =
  % V2(fid, :), where fid is all the free vertex indices. Since minFunc
  % only accepts vector variables, if V3 is a vector, we will reshape it
  % into a 3-column matrix.
  
  % cid is a column vector with all the indices of fixed variables. The
  % caller should guarantee V2(cid, :) equals to the constrained vertices
  % C.

  % The return value:
  % - arap is the arap energy, it is a scalar;
  % - darap is the column vector form of the vnum x 3 derivatives.
  %   Unfortunately we have to follow the convention from minFunc to
  %   columnize the vector before returning it.

  % Get the number of vertices.
  vnum = size(V, 1);
  
  % Build free variable index.
  fid = setdiff(1 : vnum, cid);
  
  % Reshape V3 if necessary.
  if size(V3, 2) == 1
    rows = size(V3, 1);
    V3 = reshape(V3, rows / 3, 3);
  end

  % Build V2.
  V2 = V;
  V2(fid, :) = V3;
  V2(cid, :) = C;

  % Compute E first.
  E = compe(V, V2, N, W);
  
  % Given E, a (vnum x 3) x 3 matrix, we compute R, P and S. The
  % relationship between them: for each 3 x 3 block in E:
  % [P, S] = polar(E);
  % R = P';
  % So P is the rotation matrix, S is the symmetric semi-positive definite
  % matrix. Both of P and S come from polar decomposition.
  % R is the rotation matrix in arap, which is the transpose of P.
  P = zeros(3, 3, vnum);
  S = zeros(3, 3, vnum);
  R = zeros(3, 3, vnum);
  for v = 1 : vnum
    base = (v - 1) * 3;
    [P(:, :, v), S(:, :, v)] = polar(E(base + 1 : base + 3, :));
    R(:, :, v) = P(:, :, v)';
  end
  
  % Compute arap.
  % reshape here is a little tricky. What we want for the 3rd argument is a
  % (vnum x 3) x 3 matrix, each 3 x 3 block containing the rotation matrix.
  % If we reshape P into a 3 x (vnum x 3) matrix, and P(:, i, j) will
  % become the ((j - 1) * 3 + i)-th column in the reshaped matrix. So
  % transposing it turns out to give us the correct layout for R.
  arap = comparap(V, V2, reshape(P, 3, vnum * 3)', N, W);
  
  % Now compute darap. The dimension for darap should be a vnum x 3 matrix,
  % with each entry corresponding to the variable in V2.
  darap = zeros(vnum, 3);
  
  % Before that, let's first initialize a flag to indicate which vertex is
  % free variable, and which is fixed constraints.
  isfree = ones(vnum, 1);
  isfree(cid) = 0;
  
  % Now we loop over all the edges in arap, and add each edge's
  % contributions to darap.
  for i = 1 : vnum
    % Get all the incident vertices from N.
    [~, J, ~] = find(N(i, :));
    % We first compute dEi, which is Ei's derivatives w.r.t V2. For each
    % variable, dEi is a 3 x 3 matrix. We reshape it as a 9 x 1 column
    % vector. For each vertex, we have a 9 x 3 matrix, so dEi should be a 9
    % x 3 x vnum matrix.
    dEi = zeros(9, 3, vnum);
    % Loop over all the incident vertex.
    for j = J
      % For edge(i, j), we add w * (voi - voj) * (vni - vnj)' to Ei.
      w = W(i, j);
      voi = V(i, :)';
      voj = V(j, :)';
      if isfree(i)
        dEiix = w * (voi - voj) * [1 0 0];
        dEiiy = w * (voi - voj) * [0 1 0];
        dEiiz = w * (voi - voj) * [0 0 1];
        dEi(:, :, i) = dEi(:, :, i) + [dEiix(:), dEiiy(:), dEiiz(:)];
      end

      % Compute dEijx, dEijy and dEijz.
      if isfree(j)
        dEijx = w * (voi - voj) * [-1 0 0];
        dEijy = w * (voi - voj) * [0 -1 0];
        dEijz = w * (voi - voj) * [0 0 -1];
        dEi(:, :, j) = dEi(:, :, j) + [dEijx(:), dEijy(:), dEijz(:)];
      end
    end
    
    % Get R, S, P for vertex i.
    r = R(:, :, i);
    s = S(:, :, i);
    p = P(:, :, i);
    for j = J
      w = W(i, j);
      voi = V(i, :)';
      voj = V(j, :)';
      vni = V2(i, :)';
      vnj = V2(j, :)';
      
      % Now we are focusing on the edge w * \|(vni - vnj) - r * (voi -
      % voj)\| ^ 2. 
      
      % For any variable x, the derivative can be rewritten as:
      % 2 * w * ((vni - vnj) - r * (voi - voj))' * d((vni - vnj) - 
      % r * (voi - voj)) / dx. So let's first cache it.
      cachev = 2 * w * ((vni - vnj) - r * (voi - voj))';
      
      % The next problem is how to compute d((vni - vnj) - r * (voi - voj))
      % / dx w.r.t. all the vnum x 3 variables. Note that for each variable
      % we have a 3 x 1 gradient, so we have to allocate 3 x 3 x vnum
      % space for it, with each column representing a gradient, the second
      % 3 representing the 3 elements in a single vertex.
      d = zeros(3, 3, vnum);
      % First let's consider the derivates of (vni - vnj).
      if isfree(i)
        d(:, :, i) = eye(3);
      end
      if isfree(j)
        d(:, :, j) = -eye(3);
      end
      
      % Second, we have to compute d(-r * (voi - voj)) / dx, which is
      % simply -(dr / dx) * (voi - voj).
      for k = J
        if ~isfree(k)
          continue;
        end
        for l = 1 : 3
          dei = reshape(dEi(:, l, k), 3, 3);
          dr = dpolar(dei, p, s)';
          d(:, l, k) = d(:, l, k) - dr * (voi - voj);
        end
      end
      
      % Now we have d, we can compute dv, which is this edge's derivatives
      % w.r.t V2.
      dv = reshape(cachev * reshape(d, 3, 3 * vnum), 3, vnum)';
      
      % After all is done, add dv to darap.
      darap = darap + dv;
    end
  end
  % But wait, we have to remove fixed constraints from darap.
  darap(cid, :) = [];
  darap = darap(:);
end

