function [ cots ] = tricot( T )
  % Tao Du
  % Nov 8, 2014
  
  % This function takes a triangle, returns the three cot angles.
  % The input T is a 3 x 3 matrix, each row representing one vertex.
  % The output cots is a 3 x 1 vector, each element is the corresponding 
  % cot angle of the vertex in T.
  
  % Extract the three vertices.
  A = T(1, :);
  B = T(2, :);
  C = T(3, :);
  
  % Extract the three edges.
  AB = B - A;
  BC = C - B;
  CA = A - C;
  BA = -AB;
  CB = -BC;
  AC = -CA;
  
  % Extract the length of three edges.
  a = norm(BC);
  b = norm(CA);
  c = norm(AB);

  % Allocate space for sin and cos.
  sins = zeros(3, 1);
  coss = zeros(3, 1);
  
  % Compute cos angles.
  coss(1) = AB * AC' / (c * b);
  coss(2) = BA * BC' / (c * a);
  coss(3) = CA * CB' / (b * a);
  
  % Compute sin angles.
  sins(1) = norm(cross(AB, AC)) / (c * b);
  sins(2) = norm(cross(BA, BC)) / (c * a);
  sins(3) = norm(cross(CA, CB)) / (b * a);
  
  % Compute cots.
  cots = coss ./ sins;
end