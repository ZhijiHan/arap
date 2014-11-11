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
  % V2(fid, :), where fid is all the free vertex indices.
  
  % cid is a column vector with all the indices of fixed variables. The
  % caller should guarantee V2(cid, :) equals to the constrained vertices
  % C.
  
  % We decide to loop over all the edges to manually compute arap and darap
  % for efficiency.
  
  % Get the number of vertices.
  vnum = size(V, 1);
  
  % Build free variable index.
  fid = setdiff(1 : vnum, cid);
  
  % Build V2.
  V2 = V;
  V2(fid, :) = V3;
  V2(cid, :) = C;
  
  % Get edges.
  [I, J, ~] = find(N);
  
  % Get the number of edges.
  enum = length(I);
  
  % E is a 3 x 3 x vnum temporary matrix for computing R. We initialize
  % R, P and S at the same time. The relationship between them: for each 3
  % x 3 block:
  % [P, S] = polar(E);
  % R = P';
  % So P is the rotation matrix, S is the symmetric semi-positive definite
  % matrix. Both of P and S come from polar decomposition.
  % R is the rotation matrix in arap, which is the transpose of P.
  E = zeros(3, 3, vnum);
  P = zeros(3, 3, vnum);
  S = zeros(3, 3, vnum);
  R = zeros(3, 3, vnum);
  
  % Discussion about dE: the trickest part of this function is to compute
  % E's derivatives w.r.t V2. Let's do it step by step.
  
  % First, let's fix on a single Ei = E(:, :, i) and a single V2j = V2(j,
  % :). dEi / dV2j should be a 3 x 3 x 3 matrix. Let's stack it as a 9 x 3
  % matrix, where each row k representes Ei(k)'s derivatives w.r.t. V2j.
  
  % Second, what if we have another V2j? We need another 9 x 3 matrix to
  % store it. In principle, we need a 9 x 3 x vnum matrix to hold all the
  % derivatives w.r.t all the V2.
  
  % But this is just for a single Ei! We need a 9 x 3 x vnum matrix for
  % each Ei, so finally we need to fill in a 9 x 3 x vnum x vnum matrix!
  % This is the dimension of dE:
  
  % dimension         9             3           vnum            vnum
  % explanation elements in Ei elements in V2j elements in V2 elemenst in E
  dE = zeros(9, 3, vnum, vnum);
  
  % One caveat is that we need to skip those fixed vertices. We initialize
  % a flag to determine whether a point is free or not.
  isfree = ones(vnum, 1);
  isfree(cid) = 0;
  
  % Now loop over all the edges to compute E and dE.
  for e = 1 : enum
    i = I(e);
    j = J(e);
    w = W(i, j);
    voi = V(i, :)';
    voj = V(j, :)';
    vni = V2(i, :)';
    vnj = V2(j, :)';
    
    % Update E.
    E(:, :, i) = E(:, :, i) + w * (voi - voj) * (vni - vnj)';
    
    % Update dE. Explicitly, this is updating Ei because we only update Ei.
    dEi = zeros(9, 3, vnum);
    
    % More specifically, we only need to update dEi(:, :, i) and dEi(:, :,
    % j). So dEiix is a 3 x 3 matrix, each element is Ei's derivatives w.r.t
    % V2(i, 1). dEiiy and dEiiz work similarly.
    if isfree(i)
      dEiix = w * (voi - voj) * [1 0 0];
      dEiiy = w * (voi - voj) * [0 1 0];
      dEiiz = w * (voi - voj) * [0 0 1];
      dEi(:, :, i) = [dEiix(:), dEiiy(:), dEiiz(:)];
    end
    
    % Compute dEijx, dEijy and dEijz.
    if isfree(j)
      dEijx = w * (voi - voj) * [-1 0 0];
      dEijy = w * (voi - voj) * [0 -1 0];
      dEijz = w * (voi - voj) * [0 0 -1];
      dEi(:, :, j) = [dEijx(:), dEijy(:), dEijz(:)];
    end
    
    dE(:, :, :, i) = dE(:, :, :, i) + dEi;
  end
  
  % Now given E and dE, loop over all the edges again to compute arap and
  % darap.
  
  % What is the initial value for arap?
  arap = 0.0;
  
  % Discussion about darap:
  % darap should be a vnum x 3 matrix, each element corresponds to a
  % variable in V2.
  
  % What is the initial value for darap?
  darap = zeros(vnum, 3);
  
  % Compute P, S and R first.
  for v = 1 : vnum
    [P(:, :, v), S(:, :, v)] = polar(E(:, :, v));
    R(:, :, v) = P(:, :, v)';
  end
  
  % Now we have R, let's compute arap and darap.
  for e = 1 : enum
    i = I(e);
    j = J(e);
    w = W(i, j);
    voi = V(i, :)';
    voj = V(j, :)';
    vni = V2(i, :)';
    vnj = V2(j, :)';
    r = R(:, :, i);
    
    % Remember we are now working on w * \|(vni - vnj) - r * (voi - voj)\|^2
    arap = arap + w * norm((vni - vnj) - r * (voi - voj)) ^ 2;
    
    % Now we update darap. We first want to figure out the derivatives of
    % the vector (vni - vnj) - r * (voi - voj) w.r.t. V2. This should be a
    % 3 x 3 x vnum matrix. For each V2k, there is a 3 x 3 matrix, each row
    % represents that vector's derivatives w.r.t V2k.
    dv = zeros(3, 3, vnum);
    
    % Now what should we fill in? Let's first consider (vni - vnj), which
    % should be easy: we only need to fill dv(:, :, i) and dv(:, :, j).
    if isfree(i)
      dv(:, :, i) = eye(3);
    end
    if isfree(j)
      dv(:, :, j) = -eye(3);
    end
    
    % The hard part is to figure out the influence from r. By definition
    % dE(:, :, :, i) contains all the information about r's derivatives
    % w.r.t. V2. Let's get E, P, S first:
    p = r';
    s = S(:, :, i);
    % Then Let's loop over all the V2.
    for v = 1 : vnum
      if ~isfree(v)
        continue;
      end
      % Then let's loop over all the elements in a specific V2v.
      for k = 1 : 3
        % Now what is r's derivatives w.r.t V2(v, k)?
        de = dE(:, k, v, i);
        de = reshape(de, 3, 3);
        % Now we compute dp / dV2(v, k).
        dp = dpolar(de, p, s);
        % Since we know r = p', dr / dV2(v, k) is dp's transpose.
        dr = dp';
        % Now d(-r * (voi - voj)) / dV2(v, k) is a 3 x 1 column vector.
        % This should be added into dv(:, k, v)
        dv(:, k, v) = dv(:, k, v) - dr * (voi - voj);
      end
    end
    
    % Finally, given dv, we use dv to update darap. Let's flat dv into a 3
    % x (3 x vnum) matrix, so each column is the derivative w.r.t a
    % specific V2j(k).
    dv = reshape(dv, 3, 3 * vnum);
    
    % Then we can compute the derivatives of
    % w * \|(vni - vnj) - r * (voi - voj)\| ^ 2
    % into a flat matrix.
    dflat = w * 2 * ((vni - vnj) - r * (voi - voj))' * dv;
    
    % Now dflat is a 1 x (3 x vnum) row vector, we need to reshape it into
    % a vnum x 3 matrix.
    dflat = reshape(dflat, 3, vnum);
    dflat = dflat';
    
    % Finally update darap.
    darap = darap + dflat;
  end
  
  % But wait, we need to eliminate those rows corresponding to fixed
  % vertices.
  darap(cid, :) = [];
end

