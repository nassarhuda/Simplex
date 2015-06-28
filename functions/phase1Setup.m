function [iterations,B_ids, z,x_feas] = phase1Setup(A,b)
%% This function sets up the input for simplex phase 1
% want to create E only when b<0
% minimize f'x
% st Ax <= b
% and x >= lo
bnegative = numel(find(b<0));
% phase 1 first:
% approach1, create A with an identity matrix:
[m, n] = size(A);
AI = [A,eye(m)]; % this AI includes the slack variables too
[~, nAI] = size(AI);
% check feasible sols:
% min e^T z
% subject to Ax + Ez = b, (x,z)>= 0
% where z belongs to R^m
% e = ones^T
% E is a diagonal matrix where:
% E(j,j) = +1 if b_j >= 0, -1 otherwise
% the point: (x,z) = (0,abs(b)) is feasible
% solve the system with this being the first feasible point:

% code:
% create E
E = zeros(m,bnegative);
j = 1;
B_basis = zeros(m,1);
for i=1:m
    if b(i) < 0
        E(i,j) = -1;
        B_basis(i) = nAI + j;
        j = j+1;
    else
        B_basis(i) = n + i;
    end
end

% solve the minimization problem:
chat = [zeros(nAI,1);ones(bnegative,1)];
Ahat = [AI,E];

% set up B_basis


all_ids = 1:size(Ahat,2);
N_basis = setdiff(all_ids,B_basis);

% limit is for knowing that the entering index is not bigger than it
limit = nAI;
[iterations, z, xsol, B_ids] = SimplexPhase1(Ahat,b,chat,B_basis,N_basis,limit);
x_feas = xsol(1:n);
end
