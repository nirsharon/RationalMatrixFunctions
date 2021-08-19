% figure_oscilatory_function
%
% cos(9*x)+sin(11*x) over [-1,1]
% (7,7) approximation
%
% NS, July21

clear
close all
to_save = 1;           % whether saving or not

% basic parameters
a = -1;
b = 1;
n = 7;
m = 7;
l = 1;
u = 50; 

% the function
fun = @(x) cos(9*x)+sin(11*x);

% the run
run_comparison(fun, a, b, n, m, l, u, to_save, 'oscilatory')
