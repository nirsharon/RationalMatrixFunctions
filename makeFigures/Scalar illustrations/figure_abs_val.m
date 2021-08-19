% figure_abs_val
%
% |x-.1| over [-1/2,1/2]
% (6,6)
% NS, July21

clear
close all
to_save = 1;           % whether saving or not

% basic parameters
a = -.5; %1;
b = .5;  %1;
n = 6;
m = 6;
l = 1;
u = 100;

% the function
fun = @(x) abs(x-.1);

% the run
run_comparison(fun, a, b, n, m, l, u, to_save, 'abs_val')
