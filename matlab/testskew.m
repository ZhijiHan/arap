% Tao Du
% Nov 8, 2014

% Test skew, skewmatrix and skewvector.
clc; clear all;

% Test skewmatrix.
w = rand(3, 1);
W = skewmatrix(w);

% Randomly choose x.
x = rand(3, 1);

% Verify Wx = cross(w, x).
disp('Wx = ');
disp(W * x);
disp('cross(w, x) = ');
disp(cross(w, x));

% Test skewvector.
w2 = skewvector(W);
disp('w = ');
disp(w);
disp('w2 = ');
disp(w2);

% Test skew.
A = rand(3, 3);
w = skew(A);
W = skewmatrix(w);

% Compare W with (A - A') / 2.
disp('W = ');
disp(W);
disp('(A - A^T) / 2');
disp((A - A') / 2);

clear all;