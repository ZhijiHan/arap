function [ d ] = dvi( w, vni, vnj, r, voi, voj )
  % Tao Du
  % Nov 8, 2014
  
  % Compute the derivatives of w * \|(vni - vnj) - r * (voi - voj)\|^2
  % w.r.t. vni. The return value is a 3 x 1 column vector, with vni(1), 
  % vni(2) and vni(3)'s derivatives.
  
  % Enforce all the vectors are column vectors.
  vni = vni(:);
  vnj = vnj(:);
  voi = voi(:);
  voj = voj(:);
  
  % Compute the derivatives.
  d = 2 * w * (vni - vnj - r * (voi - voj));
end

