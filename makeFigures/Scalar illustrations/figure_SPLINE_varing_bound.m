% figure_SPLINE_varing_bound
%
% same cubic spline function, fixed rational degrees, varied bound
%
% NS, July21

clear
close all
to_save = 0; %1;           % whether saving or not

% basic parameters
a = 0;
b = 3;
n = 4;
m = 5;
l = 1;

% the spline function
fun= @(x) (-x.^3 + 6*(x.^2)-6.*x+2).*((x>=0)&(x<1)) + (x.^3).*(x>=1);

% the differen cases
u1 = 2;
u2 = 4;
u3 = 8;

% run
run_comparison(fun, a, b, n, m, l, u1, to_save, 'spline1')
run_comparison(fun, a, b, n, m, l, u2, to_save, 'spline2')
run_comparison(fun, a, b, n, m, l, u3, to_save, 'spline3')
