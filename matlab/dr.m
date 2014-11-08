function [ d ] = dr( w, vni, vnj, r, voi, voj )
  % Tao Du
  % Nov 8, 2014
  
  % Compute the derivatives of w * \|(vni - vnj) - r * (voi - voj)\|^2
  % w.r.t. r. The return value d is a 3 x 3 matrix, each element in d is
  % the derivatives of the corresponding element in r.
  
  % Enforce all the vectors are column vectors.
  vni = vni(:);
  vnj = vnj(:);
  voi = voi(:);
  voj = voj(:);
  
  % Compute the derivatives.
  d = -2 * w * (vni - vnj - r * (voi - voj)) * (voi - voj)';
end

