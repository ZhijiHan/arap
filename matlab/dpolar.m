function [ dR ] = dpolar( dA, R, S )
  % Tao Du
  % Nov 8, 2014
  
  % Given A = R * S, computes R's derivatives. dA is A's derivatives.
  % Reference paper: http://run.usc.edu/substructuring/BarbicZhao-SIGGRAPH2011.pdf
  
  % Compute G.
  G = (trace(S) * eye(3) - S) * R';
  
  % Compute w.
  w = G \ (2 * skew(R' * dA));
  
  % Compute dR.
  dR = skewmatrix(w) * R;
end

