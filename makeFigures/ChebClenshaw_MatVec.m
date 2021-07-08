function [ y ] = ChebClenshaw_MatVec(A, coef, a, b, v)
    n = numel(coef);
    dim = size(A, 1);
    I = sparse(eye(dim));
    t = 2 * (A - a*I) ./ (b - a) - A; % t = 2*x;%   [-1,1]
    d = zeros(dim, 1);
    dd = zeros(dim, 1);
    for i = n : - 1 : 2
        sv = d;
        %d = 2 * t * d - dd + coef(i) * I * v;
        d = (t * d) - dd + coef(i) * v;
        dd = sv;
    end
    y = A * d - dd + 0.5 * coef(1) * v;
end
