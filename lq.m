% Defining a function such L,Q will be returned  
function [L, Q] = lq(matrixA)
% Setting n for when we do our function
[~, n] = size(matrixA);
% Set Q and L
Q = eye(n);
L=matrixA;
% Defining our loop such that we are looking at rows one by one and
% applying the transformation to each row
for k=1:size(L,1)
    ckt=L(k,k:end);
    L(:,k:end)=rhouse(ckt,(L(:,k:end)));
    Q(:,k:end)=rhouse(ckt,(Q(:,k:end)));
end
L;
Q=Q';
% Gives the option of observing the points plotted in matrix L
spy(abs(L)>1e-14);
end

% Local function to apply Householder reflection which has been adapted
% from Lab2 solutions
function w=rhouse(c,w)
% Householder matrix-vector product
% Input: column vectors c,w
% Output: the vector w is overwritten by Q(u(c'))*w
% c has been transposed to target the rows of the matrix
u=housevec(c');
% Adjusted line to fit the targeting of the rows  
w=w-2*w*(u*u');
end

function u=housevec(c)
% This file generates the Householder unit vector u(c)
% such that Q(u(c))c=norm(c)e1.
e1=eye(length(c),1);
u=c-norm(c)*e1;
u=u/norm(u);
end