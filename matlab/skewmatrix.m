function [ W ] = skewmatrix( w )
  % Tao Du
  % Nov 8, 2014
  
  % Given w, a 3 x 1 column vector, computes a 3 x 3 skew matrix W, such
  % that Wx = cross(w, x) for all the 3 x 1 column vector x.
  
  % Initialize W.
  W = zeros(3, 3);
  
  % Compute each column in W.
  W(:, 1) = cross(w, [1; 0; 0]);
  W(:, 2) = cross(w, [0; 1; 0]);
  W(:, 3) = cross(w, [0; 0; 1]);
end

