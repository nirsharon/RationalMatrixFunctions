function fv=cheb_eval_mat_vec(coef, A, m, a, b, v)
% Chebyshev evaluation: The Chebyshev polynomial
%       \sum_{k=0}^{m-1} c_{k}T_{k}(A) - c_{0}/2
%
% Translated from Numerical Recipes, Third edition, Section 5.8, pp. 237.
%
% NS, Dec 19 (a joint work with Y.Shkolnisky) + August 21

if nargin<4
    a=-1;
    b=1;
end
d=zeros(size(v));
dd=zeros(size(v));

% if issparse(A)
%     if max(eigs(A))>b && min(eigs(A))<a
%         error('The matrix does not fit the approximation');
%     end
% else
%     if norm(A)>1
%         error('The matrix does not fit the approximation');
%     end
% end


if m>numel(coef) % Lowest order of apprixmation is one, which corrponds to the constant approximation
    error('Approximation order is too high for the precomputed coefficients');
end

if m<1
    error('Approximation order must be greater than 1');
end

if (a~=-1)&&(b~=1)
    n = size(A,1);
    newA=(2*A+(-a-b)*eye(n))*((1/(b-a))*eye(n));
else
    newA = A;
end

y2=2*newA;
for j=m-1:-1:1 % Clenshaw's recurrence.
    sv=d;
    d=y2*d-dd+coef(j+1)*v;
    dd=sv;
end
% fv = newA*(d*v)-dd*v+0.5*coef(1)*v;
fv = newA*(d)-dd+0.5*coef(1)*v;
end
