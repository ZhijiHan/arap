function [ d ] = dvj( w, vni, vnj, r, voi, voj )
  % Tao Du
  % Nov 8, 2014
  
  % Compute the derivatives of w * \|(vni - vnj) - r * (voi - voj)\|^2
  % w.r.t. vnj. The return value is a 3 x 1 column vector, with vnj(1), 
  % vnj(2) and vnj(3)'s derivatives.
  
  % Enforce all the vectors are column vectors.
  vni = vni(:);
  vnj = vnj(:);
  voi = voi(:);
  voj = voj(:);
  
  % Compute the derivatives.
  d = -2 * w * (vni - vnj - r * (voi - voj));
end

