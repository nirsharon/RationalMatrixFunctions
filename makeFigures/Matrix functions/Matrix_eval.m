% Matrix rational approximation evaluation
% Based on chebeval_matrix function
% Nir Sharon, Elior Kalfon, March 2020.
function Mrat=Matrix_eval(p,q,A,a,b,condtest)

% Evaluate matrix using optimization method 
% Input:
%  p - numerator coefficients
%  q - denominator coefficients
%  A - matrix to be evaluated
%  a - region lower bound
%  b - region upper bound
%  condtest - condition value test
% Ouput:
%   Mrat - the evaluated matrix 

if nargin==3
    a=-1;
    b=1;
end
if nargin<6
    condtest=0;
end
if condtest==0
    %Cheb evaluation of numerator and denominator
    Tp=chebeval_matrix(p,A,length(p),a,b);
    Tq=chebeval_matrix(q,A,length(q),a,b);
    Mrat = Tp/Tq;
    %cond(Tq);
else
    %Cheb evaluation of denominator
    q(1) = 2*q(1);
    Tq=chebeval_matrix(q,A,length(q),a,b);
    Mrat = inv(Tq);
end