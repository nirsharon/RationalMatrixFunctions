function fv=cheb_eval_vector_AxFunc(coef, AxFunc ,m , v)
% Chebyshev evaluation: The Chebyshev polynomial
%       \sum_{k=0}^{m-1} c_{k}T_{A}v - c_{0}/2*v
% We assume A is not given explicitly, but can be evaluated as 
% A(x) = AxFunc(x), where x is a vector
%
% Translated from Numerrical Recipes, Third edition, Section 5.8, pp. 237.
%
% July 2016. Revised 2021 NS
n  =  size(v,1);
l  =  size(v,2);
d  = zeros(n,l);
dd = zeros(n,l);

%y2=2*A; 

if m>numel(coef) % Lowest order of apprixmation is one, which corrponds to the constant approximation
    error('Approximation order is too high for the precomputed C');
end

if m<1
    error('Approximation order must be greater than 1');
end

for j=m-1:-1:1 % Clenshaw recurrence.
    sv=d;
    d=2*AxFunc(d)-dd+coef(j+1)*v;
    dd=sv;
end
fv= AxFunc(d)-dd+0.5*coef(1)*v;
end
