% Tao Du
% Nov 8, 2014

% Test dpolar.
clc; clear all;

% Randomly generate a matrix A.
A = rand(3, 3);

% Display the determinate of A.
disp('det(A) = ');
disp(det(A));

% A can be treated as A(t = 0) = A + t, so dA = ones(3, 3).
dA = ones(3, 3);

% Compute polar decomposition.
[R, S] = polar(A);

% Now add a little perturbation.
epsilon = 1e-4;
A2 = A + epsilon;

% Compute the new polar decomposition.
[R2, ~] = polar(A2);

% Compute the numerical derivatives.
dRn = (R2 - R) / epsilon;

% Compute the analytic derivatives.
dRa = dpolar(dA, R, S);

% Compare the results.
disp('Numerical derivatives = ');
disp(dRn);
disp('Analytic derivatives = ');
disp(dRa);
disp('Difference = ');
disp(dRn - dRa);

clear all;
