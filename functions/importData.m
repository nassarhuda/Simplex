function[A,b,c] = importData(fileName)
% imports data from the mat file 
load(fileName);
A = Problem.A;
b = Problem.b;
% c = Problem.aux.c;
[m,n]=size(A);
c = ones(n,1); % handout requirement
%lo = Problem.aux.lo;
%hi = Problem.aux.hi;










