% figure_cusp
%
% |x|^2/3 over [-1,2]
% (6,4)
% NS, July21

clear
close all
to_save = 0;           % whether saving or not

% basic parameters
a = -1;
b = 2;
n = 4;
m = 4;
l = 1;
u = 100;
    

% the function
fun = @(x) abs(x).^(2/3);

% the run
run_comparison(fun, a, b, n, m, l, u, to_save, 'cusp')
