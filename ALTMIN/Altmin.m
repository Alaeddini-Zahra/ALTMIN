
clear;
clc;
% ///////////// Parameter Setting /////////////////
m = 200;
n = 200;
r = 3; % Rank of matrix
C = 4; % Specify number of samples: C times of dimension.

NumRange = 10; % Range of random numbers that generate the orignal matrix
maxIter = 30; % Maximum number of iterations for AltMin method
tol = 0.0001; % Tolerance of difference below which AltMin will stop iteration
show = true; % Iteration log will show if this is true



% ////////////// Running Area ///////////////////
NumOfScale = m * n;
NumOfDimension = m*r + n*r - r^2;
p = C * NumOfDimension / NumOfScale; % Sampling probability

a = randi(NumRange,m,r);
b = randi(NumRange,r,n);
A = a*b;
MASK = zeros(m,n);
MASK = RandomPick(MASK,p);
MASK_ = ~MASK;
NumOfSample = sum(MASK(:));

AA = A .* MASK;

tic;
[Acomplete,difference] = AltMin3(AA,MASK,r,maxIter,tol,show);
toc;
diffALL = norm(MASK_.*(A - Acomplete),'fro')/norm(MASK_.*A,'fro'); % for unmasked entries
