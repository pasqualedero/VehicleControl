clear; close all; clc;

%% HOW TO SOLVE LMI? 

% ------------------------- YALMIP --------------------------------------



% ------------------------- LMI TOOLBOX ----------------------------------

%setlmis([]); % beginning of LMI system





%% Exercise 1
% Solve the following optimization problem, with Yalmip + solver (SeDuMi or LMIlab)

A = [-1 5; 0 -2];
B = [0 1]';
Cm = [1 1];
D = -27;

[n,m] = size(A);
[p,q] = size(D);

epsilon = 10e-8;

% LMI def.
X = sdpvar(n);      % symmetric n,n matrix
gamma = sdpvar(1);  % scalar

% Constraints
C1 = X >= epsilon;
C2 = gamma >= epsilon;
C3 = ([A'*X+X*A X*B     Cm' ;
        (X*B)' -gamma   D' ;
        Cm       D       -gamma] <= epsilon);

C = C1+C2+C3;

opts = sdpsettings;
opts.solver = 'sedumi'; % lmilab

solution = optimize(C,gamma,opts);
Xsol = double(X);
gammaSol = double(gamma);

matlab = normhinf(A,B,Cm,D); % matlab command for H_inf norm

