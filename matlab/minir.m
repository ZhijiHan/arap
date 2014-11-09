function [ R ] = minir( V, V2, N, W )
  % Tao Du
  % Nov 8, 2014
  
  % Given vertices, new vertices, neighborhood and weight, minimize
  % rotation matrices by SVD projection.
  
  % The return value R is a (#vertices x 3) x 3 matrix.
  
  % Get # vertices.
  vnum = size(V, 1);
  
  % Compute E.
  E = compe(V, V2, N, W);
  
  % Preallocate space for rotation matrix.
  R = zeros(3 * vnum, 3);
  
  % Use SVD to solve rotation matrix.
  for i = 1 : vnum
    % Compute the base index for vertex i.
    base = (i - 1) * 3;
    
    % Compute r by polar decomposition.
    [r, ~] = polar(E(base + 1 : base + 3, :));
    r = r';
    
    % Write r back to R.
    R(base + 1 : base + 3, :) = r;
  end
end

