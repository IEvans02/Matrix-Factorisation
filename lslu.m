function [x,L_E,U_E] = lslu(A,b)
%function to implement the LU Factorisation method for solving 
%the Least Squares problem for m>n
%Input: matrix A and vector b
%Output: solutions to linear system x and economical matrices L_E and U_E
m=size(A,1);
n=size(A,2);
A=ludecompgen(A);
U = triu(A);
L = (tril(A,-1)*eye(n,m))+eye(m);

%Setting L_E and U_E
L_E = L*eye(m,n);
U_E = eye(n,m)*U;
 
B = L_E'*L_E;
%factorise B=R'*R using Cholesky to find R
R=mychol(B);
%Then R'(R*ustar)=LE'*b
%Set C:=LE'*b
C=L_E'*b;
%Solve for z:=R*ustar
z = forwardsub(R',C);
%Solve for ystar: R.*ystar = z
%y = inv(R)*z;
y = backwardsub(R,z);
%Solve for xstar: UE.*xstar = ystar
x = backwardsub(U_E,y);
end%lslu


function A=ludecompgen(A)
%LU factorisation of a dense matrix A such that m>n.
%Input: dense matrix A with m>n.
%Output: matrix A overwritten by economical matrices LE and UE.
n=size(A,2);
m=size(A,1);
p=min(m,n);
if m==n 
    p=p-1;
end    
for k=1:p
    for j=k+1:m
        if A(k,k)==0, error('Zero pivot: LU factorisation does not exist.');end
        A(j,k)=A(j,k)/A(k,k);
        for i=k+1:n
            A(j,i)=A(j,i)-A(j,k)*A(k,i);
        end
    end
end
end%ludecompgen


function A=mgauss(c,A)
%This function computes the Gauss transformation L(c) of A
%Input:  column vector c and matrix A
%Output: the matrix A is overwritten by the matrix L(c)*A

%Test for existence of L(c):
if c(1)==0
    error('Gauss transformation not defined')
end
A(2:end,:)=A(2:end,:)-c(2:end,1)*A(1,:)/c(1);
end%mgauss


function x=backwardsub(U,b)
%function x=backwardsub(U,b)
%BACKWARDSUB This function file computes the solution of Ux=b
%using the backward substitution algorithm (see Examples sheet 1, Q11).
%Input: upper triangular matrix U, column vector b
%Output: solution x of linear system Ux=b.
n=length(b);
x=zeros(n,1);
x(n)=b(n)/U(n,n);
for i = n-1:-1:1
    x(i)=b(i)-U(i,i+1:n)*x(i+1:n);
    x(i)=x(i)/U(i,i);
end
end%backwardsub


function A = mychol(A)
%factorise for R such that A=R'*R
%Input: matrix A such that A = R'*R for some R we seek
%Output: factorised matrix R
for k=1:length(A)
    ctk=A(k:end,k);
    A(k:end,:)=mgauss(ctk,A(k:end,:));
    A(k,:) = A(k,:)/sqrt(A(k,k));
end    
end%mychol


function x=forwardsub(L,b)

%function x=forwardsub(L,b)
%FORWARDSUB This function file computes the solution of Lx=b
%using the forward substitution algorithm.
%Input: square lower triangular matrix L, column vector b
%Output: solution x of linear system Lx=b.
n=length(b);
x=zeros(n,1);
x(1)=b(1)/L(1,1);
for i=2:n
    x(i)=b(i)-L(i,1:i-1)*x(1:i-1);
    x(i)=x(i)/L(i,i);
end

end%forwardsub