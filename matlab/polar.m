function [ R, S ] = polar( A )
  % Tao Du
  % Nov 8, 2014
  
  % Given a matrix A, returns polar decompositions R, S such that R * S = A
  % and R is orthogonal matrix with det(R) = 1, and S is a symmetric
  % matrix.
  % Strictly speaking, this is not the correct polar decomposition because
  % the actual polar decomposition requires S be a symmetric positive
  % semi-definite matrix, and det(R) can be 1. However, in our example we
  % decide to use the definition above because we require R be a rotation
  % matrix.
  
  % Do SVD for A.
  [u, sig, v] = svd(A);
  
  % Now we know A = u * sig * v', let's compute R = u * v' first. If
  % det(R) < 0, we flip the last column in u, and change the sign of the
  % last singular value in sig.
  R = u * v';
  if det(R) < 0
    u(:, end) = -1 * u(:, end);
    sig(end, end) = -1 * sig(end, end);
    R = u * v';
  end
  
  % Note that A = u * sig * v' still holds, we can compute S as v * sig *
  % v'.
  S = v * sig * v';
end

