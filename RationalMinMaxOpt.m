function [p, q, zval] = RationalMinMaxOpt(f, n, m, pts, LB, UB, a, b, prc, vrb)
% Calculating the uniform best rational approx of type (n,m)
% via optimization with deviation precision 'eps'
%
% Input:
%   f   - the function to be approximated. A function handler
%   n,m - the rational approx parameters = maximum degree (numer.,deno.)
%   pts - discretization points
%   LB  - lower bound on the denominator (away from zero)
%   UB  - upper bound on the denominator
% 	prc - precision of the bisection (maximum deviation accuracy)
%   vrb - flag for verbose run
% Ouput:
%   p,q - the rational approx coefficients
%   z   - the maximal deviation
% Run example:
%  [p, q, z] = RationalMinMaxOpt(@(x) abs(x), 4, 4, linspace(-1,1,20), .1, 100, 1e-10, 1)
%
% Elior Kalfon, Nir Sharon, Feb 2020

% verbose run
if nargin < 10
    vrb = 0;
end
% deviation precision
if nargin < 9
    prc = 10^-15;
end

if nargin < 7
    a=-1;
    b=1;
end
% upper bound on denominator (not necessary in general)
if nargin < 6
    UB = 1000*LB;
end

% "Chebyshev" Vandermonde matrix
I = eye(n); I(1) = 2;
Tn = zeros(numel(pts),n);
for deg=0:(n-1)
    Tn(:,deg+1) = chebeval_scalars(I(deg+1,:), pts, deg+1,a,b);
end
I = eye(m); I(1) = 2;
Tm = zeros(numel(pts),m);
for deg=0:(m-1)
    Tm(:,deg+1) = chebeval_scalars(I(deg+1,:), pts, deg+1,a,b);
end

% lower bound
uL  = 0;

bma=0.5*(b-a);
bpa=0.5*(b+a);

% upper bound by polynomial interpolation (degree n+m, chebyshev pts)
numpts = n+1; 
int_pts = vec(bpa+bma*cos( pi* (2.*( numpts:-1:1) -1 ) / (2*numpts )));
uH = max(abs( f(pts) - barycentric_poly_interpolation(int_pts, f(int_pts), pts)));
while(uH-uL)>prc
    z=(uH+uL)/2;
    if(checkVal(f,z,n,m,pts,Tn,Tm,LB,UB,vrb,inf,-inf))
        uH=z;
    else
        uL=z;
    end
end

%Calculates the optimal p,q
[p,q,zval]=LpRat(f,uH,n,m,pts,Tn,Tm,LB,UB,vrb,inf,-inf);
end

function bool=checkVal(f,z,n,m,T,Tn,Tm,DLB,DUB,vrb,uH,uL)
[~,~,u,~]= LpRat(f,z,n,m,T,Tn,Tm,DLB,DUB,vrb,uH,uL);
bool     = (u<= 1e-15);
end

function [p,q,u,exitflag] = LpRat(f,z,n,m,T,Tn,Tm,DLB,DUB,vrb,uH,uL)

%cond1 and cond2 are the one by column multiplication of each site and a row of
%chebyshev matrix
f1 = f(T)+z;
f2 = f(T)-z;
cond1 = f1(:).*Tm(:,:);
cond2 = f2(:).*Tm(:,:);

% A is a matrix of constraints corresponds to the inequation Ax<b
A=[Tn -cond1 -ones(length(T),1);
    -Tn cond2 -ones(length(T),1);
    zeros(size(Tn)) -Tm zeros(length(T),1);
    zeros(size(Tn)) Tm zeros(length(T),1)];

b=[zeros(2*size(Tn,1),1);
    -DLB*ones(length(T),1);
    DUB*ones(length(T),1)];

lb=[-inf.*ones(1,n) -inf -inf.*ones(1,m-1) uL];
ub=[inf.*ones(1,n) inf  inf.*ones(1,m-1) uH];

%A dummy varibale in order to arrange in one objective function
obj=[zeros(1,n+m) 1];

if vrb
    options = optimoptions('linprog','Display','iter');
else
    options = optimoptions('linprog','Display','None','algorithm','dual-simplex');
end
if vrb
    [x,~,exitflag,output] = linprog(obj,A,b,[],[],lb,ub,options);
    output
else
    options.ConstraintTolerance = 1e-9;
    options.OptimalityTolerance = 1e-9;
    [x,~,exitflag,outm] = linprog(obj,A,b,[],[],lb,ub,options);
end
if isempty(x)
    p=0;
    q=0;
    u=inf;
else
    p=x(1:n);
    q=x(n+1:m+n);
    u=x(n+m+1);
end
end

