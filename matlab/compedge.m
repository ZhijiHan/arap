function [ e ] = compedge( w, vni, vnj, r, voi, voj )
  % Tao Du
  % Nov 8, 2014
  
  % Compute the energy for each edge.
  % e = w * \|(vni - vnj) - r * (voi - voj)\|^2.
  
  % Enforce all the vectors are column vectors.
  vni = vni(:);
  vnj = vnj(:);
  voi = voi(:);
  voj = voj(:);
  
  % Compute the edge energy.
  e = w * norm((vni - vnj) - r * (voi - voj)) ^ 2;
end

