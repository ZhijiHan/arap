function [ w ] = skewvector( W )
  % Tao Du
  % Nov 8, 2014
  
  % Given a 3 x 3 skew matrix W, computes its corresponding skew vector w.
  % This is the inverse function of skewmatrix. i.e.
  % skewvector(skewmatrix(w)) = w.
  % All we know is that Wx = cross(w, x) for all the 3 x 1 column vector x.
  % If we denote w = [w1; w2; w3], W can be represented as
  % W = [0 -w3 w2; w3 0 -w1; -w2 w1 0].
  
  % Initialize w.
  w = zeros(3, 1);
  w(1) = W(3, 2);
  w(2) = W(1, 3);
  w(3) = W(2, 1);
end

