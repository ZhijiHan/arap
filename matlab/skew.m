function [ w ] = skew( A )
  % Tao Du
  % Nov 8, 2014
  
  % Computes the unique skew-vector w so that skewmatrix(w) = (A - A') / 2.
  w = skewvector((A - A') / 2);
end

