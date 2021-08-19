% run_examples_script

% rational approximation playground 

clear
close all

run_case = 5; %[1,2,3,4,5];
to_save = 0; %1;           % whether saving or not

% case 1: spline 
if ismember(1,run_case)
    close all

    a = 0;
    b = 3;
    n = 4;
    m = 5;
    l = 1;
    u = 2;
    fun= @(x) (-x.^3 + 6*(x.^2)-6.*x+2).*((x>=0)&(x<1)) + (x.^3).*(x>=1);
    run_comparison(fun, a, b, n, m, l, u, to_save, 'spline')
end

% case 2: oscilatory function
if ismember(2, run_case)
    close all

    a = -1;
    b = 1;
    n = 7;
    m = 7;
    l = 1;
    u = 50;
    fun = @(x) cos(9*x)+sin(11*x);
    run_comparison(fun, a, b, n, m, l, u, to_save, 'oscilatory')
end

% case 3: high derivative values
if ismember(3,run_case)
    close all

    a = -1;
    b = 1;
    n = 4;
    m = 10;
    l = 1;
    u = 50;
    fun = @(x) x.^20;
    run_comparison(fun, a, b, n, m, l, u, to_save, 'high_dev')
end

%====>
% case 4: jump in derivative
if ismember(4,run_case)
    close all

    a = -.5; %1;
    b = .5;  %1;
    n = 6;
    m = 6;
    l = 1;
    u = 100;
    fun = @(x) abs(x-.1);
    run_comparison(fun, a, b, n, m, l, u, to_save, 'abs_val')
end

% case 5: cusp function
if ismember(5, run_case)
    close all
    
    a = -1;
    b = 2;
    n = 6;
    m = 4;
    l = 1;
    u = 100;
    fun = @(x) abs(x).^(2/3);
    run_comparison(fun, a, b, n, m, l, u, to_save, 'cusp')
end
