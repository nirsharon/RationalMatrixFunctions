% error_denominator_script
% case of cubic spline

clear
tic

% main parameters
to_save = 1;
degree_vals = 5:11; %4:2:10;
a = 0;
b = 3;
l = 1;
u = 100;

% the function
fun = @(x) (-x.^3 + 6*(x.^2)-6.*x+2).*((x>=0)&(x<1)) + (x.^3).*(x>=1);
f1  = chebfun(fun,[a,b]);
eps = 1e-15;

% error evaluation
ev_N = 10^3;
ev_pts = linspace(a, b, ev_N);
ev_pts = ev_pts(:);

% sampling
N   = 512; 
pts = linspace(a, b, N);
pts = pts(:);

err_remez = zeros(size(degree_vals));
den_remez = zeros(size(degree_vals));
err_opt = zeros(size(degree_vals));
den_opt = zeros(size(degree_vals));

for j=1:numel(degree_vals)
    n = degree_vals(j);
    [~, q_poly, r] = minimax(f1, n-1, n);
    err_remez(j) = max(abs(f1(ev_pts)-r(ev_pts)));
    den_remez(j) = max(abs(q_poly(ev_pts)))/min(abs(q_poly(ev_pts)));
    
    [p, q, max_dev] = RationalMinMaxOpt(fun, n, n+1, pts, l, u, a, b, eps);
    p(1) = 2*p(1); q(1) = 2*q(1);
    Tp   = chebeval_scalars(p, ev_pts, n, a, b);
    Tq   = chebeval_scalars(q, ev_pts, n+1, a, b);
    app  = Tp(:)./Tq(:);
    err_opt(j) = max(abs(app - fun(ev_pts)));
    den_opt(j) = max(abs(Tq(:)))/min(abs(Tq(:)));
    
end

%%

figure;
set(0,'defaultTextInterpreter','latex');

yyaxis left
h1 = semilogy(degree_vals,err_opt,'LineWidth',3.5);
hold on
h2 = semilogy(degree_vals,err_remez,'--','LineWidth',3);
xlabel('The degree m')
ylabel('Uniform Error')

yyaxis right
ylabel('$C_r$')
set(0,'defaultTextInterpreter','latex');

h3 = semilogy(degree_vals, den_opt,'LineWidth',3);
h4 = semilogy(degree_vals,den_remez,'--','LineWidth',3);
set(gca,'FontSize',18)

% Produce left-axis legend
leg1=legend([h1 h2],'Optimization', 'Remez', 'Location', 'best');

% Create invisible axis in the same position as the current axes
ax = gca(); % Handle to the main axes
axisCopy = axes('Position', ax.Position, 'Visible', 'off');

% Copy objects to second axes
hCopy3 = copyobj(h3, axisCopy);
hCopy4 = copyobj(h4, axisCopy);

% Replace all x values with NaN so the line doesn't appear
hCopy3.XData = nan(size(hCopy3.XData));
hCopy4.XData = nan(size(hCopy4.XData));

% Create right axis legend
legend([hCopy3, hCopy4],'Optimization','Remez', 'Location', 'SouthWest')
set(leg1,'FontSize',18);
set(gca,'FontSize',18)

if to_save
    folder_name = ['error_denom_',datestr(now,'mmmm_dd_yy')];
    mkdir(folder_name)
    cd(folder_name)
    
    name_it = 'error_denominator';
    saveas(gcf, name_it ,'fig');
    saveas(gcf, name_it,'jpg');
    print('-depsc2',name_it);
    
    save('error_denom_data');
    cd '../'
end
toc()

