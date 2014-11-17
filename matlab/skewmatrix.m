function [ W ] = skewmatrix( w )
  % Tao Du
  % Nov 8, 2014
  
  % Given w, a 3 x 1 column vector, computes a 3 x 3 skew matrix W, such
  % that Wx = cross(w, x) for all the 3 x 1 column vector x.
  
  W = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

