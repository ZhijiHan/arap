% Tao Du
% Nov 8, 2014

% Testing the computation of derivatives in dvi, dvj and dr.

% Initialize w, vni, vnj, r, voi, voj randomly.
clc; clear all;
w = rand();
vni = rand(3, 1);
vnj = rand(3, 1);
r = rand(3, 3);
voi = rand(3, 1);
voj = rand(3, 1);

% Compute the energy.
e = compedge(w, vni, vnj, r, voi, voj);

% Perturbe the element a little bit.
epsilon = 1e-4;

%% Test dvi.
% Compute the dvi numerically.
dvni = zeros(3, 1);
for i = 1 : 3
  % Add perturbations in i-th element.
  vni(i) = vni(i) + epsilon;
  
  % Compute the new energy.
  e2 = compedge(w, vni, vnj, r, voi, voj);
  
  % The numerical derivatives.
  dvni(i) = (e2 - e) / epsilon;
  
  % Reset i-th element.
  vni(i) = vni(i) - epsilon;
end

% Get the analytic solution.
dvai = dvi(w, vni, vnj, r, voi, voj);

% Display the two solutions.
disp('numerical solution:');
disp(dvni);
disp('analytic solution');
disp(dvai);

%% Test dvj.
% Compute the dvi numerically.
dvnj = zeros(3, 1);
for i = 1 : 3
  % Add perturbations in i-th element.
  vnj(i) = vnj(i) + epsilon;
  
  % Compute the new energy.
  e2 = compedge(w, vni, vnj, r, voi, voj);
  
  % The numerical derivatives.
  dvnj(i) = (e2 - e) / epsilon;
  
  % Reset i-th element.
  vnj(i) = vnj(i) - epsilon;
end

% Get the analytic solution.
dvaj = dvj(w, vni, vnj, r, voi, voj);

% Display the two solutions.
disp('numerical solution:');
disp(dvnj);
disp('analytic solution');
disp(dvaj);

%% Test dr.
% Compute dr numerically.
dnr = zeros(3, 3);
for i = 1 : 3
  for j = 1 : 3
    % Add perturbations in i-th row and j-th column.
    r(i, j) = r(i, j) + epsilon;
    
    % Compute the new energy.
    e2 = compedge(w, vni, vnj, r, voi, voj);
    
    % The numerical derivatives.
    dnr(i, j) = (e2 - e) / epsilon;
    
    % Reset the element.
    r(i, j) = r(i, j) - epsilon;
  end
end

% Get the analytic solution.
dar = dr(w, vni, vnj, r, voi, voj);

% Display the two solutions.
disp('numerical solution:');
disp(dnr);
disp('analytic solution');
disp(dar);

clear all;