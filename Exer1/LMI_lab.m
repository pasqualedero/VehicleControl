clear; close all; clc;

%% Exercise 1
% Solve the following optimization problem, with Yalmip + solver (SeDuMi or LMIlab)

A = [-1 5; 0 -2];
B = [0 1]';
C = [1 1];
D = -27;

[n,m] = size(A);
[p,q] = size(D);

epsilon = 10e-8;

setlmis([]);
X = lmivar(1, [2 1]);
gamma = lmivar(1, [1 1]);

i = 1;
lmiterm([-i 1 1 X],1,1);    % X>0
i=i+1;
lmiterm([-i 1 1 gamma],1,1);  % gamma > 0  
i=i+1;
lmiterm([i 1 1 X],1,A,'s');
lmiterm([i 1 2 X],1,B);
lmiterm([i 1 3 0],C');
lmiterm([i 2 2 gamma],-1,1);
lmiterm([i 1 2 X],1,B);
lmiterm([i 2 3 0],D');
lmiterm([i 3 3 gamma],-1,1);

energy = getlmis;
N = decnbr(energy);

for j=1:N
    c(j) = defcx(energy,j,gamma);
end

[copt, xopt] = mincx(energy, c);
