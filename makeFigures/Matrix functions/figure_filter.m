% figure_filter
%
% This script runs the matrix filtering example 
%
% NS, EK, July 21

clear
close all
to_save = 0;

%% main parameters
a = -1;
b = +1;
    
n = 10; 
m = 10; 
    
l = 1;
u = 1000;

trial_name = 'filtering_matrix_function';

%% filter function
r = .05;    % rise
c = .4;    % center
R = .2;     % Radius
fun = @(x) x.*(.5*(1 - erf((2/r)*(abs(x-c) - R ) )));


%% sampling points
N   = 500; 
pts = linspace(a, b, N);
pts = pts(:);

% error evaluation
ev_N = 1+10^3;
ev_pts = linspace(a, b, ev_N);
ev_pts = ev_pts(:);

%% our approximation
tic; 
eps = 1e-15;
[p, q, ~] = RationalMinMaxOpt(fun, n+1, m+1, pts, l, u, a, b, eps); 
time_opt = toc;

% evaluate the result
p(1) = 2*p(1);
q(1) = 2*q(1);
Tp   = chebeval_scalars(p, ev_pts, n+1, a, b);
Tq   = chebeval_scalars(q, ev_pts, m+1, a, b);
app  = Tp(:)./Tq(:);
err_opt = (app - fun(ev_pts));

%% polynomial Remez approximation 
warning('off')
f1 = chebfun(fun,[a,b]);
total_deg = n+m;
p_remez = minimax(f1, total_deg); 
warning('on')

poly_app   = p_remez(ev_pts);
err_remez  = f1(ev_pts)-poly_app;
poly_coefs = polyfit(ev_pts,p_remez(ev_pts),total_deg);

%% triple-A 
max_iters = 200;
[r1, pol, res, zer, z, f, w, errvec] = aaa(fun(pts), pts, 'degree', n, 'lawson', max_iters);
aaa_ap  = r1(ev_pts);
err_aaa = fun(ev_pts)-r1(ev_pts);

% denominator of AAA
aaa_pts = setdiff(ev_pts, z);
cauchy_mat = 1./bsxfun(@minus,aaa_pts,z.'); 
aaa_denom  = abs(cauchy_mat*w.*(prod(aaa_pts'-z)'));

% coefs
q_aaa = polyfit(aaa_pts, aaa_denom, m);
p_aaa = polyfit(aaa_pts, aaa_denom.*r1(aaa_pts), n);

%% matrix function
N = 100;            % matrix size
Q = RandOrthMat(N); % Rand orthogonal matrix
spect = vec(0.5*(b+a)+0.5*(b-a)*cos( pi* (2.*( N:-1:1) -1 ) / (2* N) ));

% evaluate in SINGLE precision
A_sing   = single(Q*diag(spect)*Q');
A   = (Q*diag(spect)*Q');
gtA = Q*diag(fun(spect))*Q';

rA    = Matrix_eval(p, q, A, a, b, 0);
aaarA_single = single(polyvalm(p_aaa, A_sing))/single(polyvalm(q_aaa, A_sing)); % Horner based
aaarA = (polyvalm(p_aaa, A))/(polyvalm(q_aaa, A)); % Horner based
polyA = polyvalm(poly_coefs, A); 

err_opt_mat  = norm(rA - gtA,'fro')/norm(gtA,'fro');
err_aaa_mat  = norm(aaarA - gtA,'fro')/norm(gtA,'fro');
err_poly_mat = norm(polyA - gtA,'fro')/norm(gtA,'fro');
% cond(polyvalm(q_aa,A))

% open a folder if we need to save
if to_save
    folder_name = [trial_name,'_',datestr(now,'mmmm_dd_yy')];
    mkdir(folder_name)
    cd(folder_name)
end

%% plot the filter
figure;
set(0,'defaultTextInterpreter','latex');
f_pts = linspace(-1,1,500);
plot(f_pts, fun(f_pts),'linewidth', 3);
hold on;
s1 = c-(R-.7*r);
e1 = c+(R-.9*r);
plot(linspace(s1, e1, 10), linspace(s1, e1, 10),'--r','linewidth', 3);
plot(f_pts, f_pts.*( heaviside(f_pts - s1) - heaviside(f_pts - e1)) )
legend({'The filter function', 'Effective zone','Ideal filter'},'Location','NorthWest','Interpreter','latex')
set(gca,'FontSize',18)

if to_save
    nameit = [trial_name,'_filter_fun'];
    saveas(gcf, nameit ,'fig');
    saveas(gcf, nameit,'jpg');
    print('-depsc2',nameit);
end

%% plot the histogram
hist_bins = 20;
bin_edges = linspace(-1,1,hist_bins); 
set(0,'defaultTextInterpreter','latex');
figure; h1 = histogram(spect,bin_edges); 
hold on; h2 = histogram(fun(spect),bin_edges);
legend('Original','Filtered')
set(gca,'FontSize',18)

if to_save
    nameit = [trial_name,'_hists'];
    saveas(gcf, nameit ,'fig');
    saveas(gcf, nameit,'jpg');
    print('-depsc2',nameit);
end

%% the function approximations
figure
set(0,'defaultTextInterpreter','latex');
plot(ev_pts, app,'linewidth', 3);
hold on;
plot(ev_pts, poly_app,'--r','linewidth',3.5);
plot(ev_pts, aaa_ap,':b','linewidth',3.5);
legend({'Optimization','Remez polynomial','AAA'},'Location','NorthWest','Interpreter','latex')
grid on
set(gca,'FontSize',18)

if to_save
    nameit = [trial_name,'_approxs'];
    saveas(gcf, nameit ,'fig');
    saveas(gcf, nameit,'jpg');
    print('-depsc2',nameit);
end

%% denominaor values
figure
semilogy(ev_pts, Tq(:),'linewidth', 3);
hold on;
semilogy(aaa_pts, aaa_denom,':b','linewidth',3.5);
legend({'Optimization','AAA'},'Location','SouthWest','Interpreter','latex')
grid on
set(gca,'FontSize',18)
c1 = max(abs(Tq(:)))/min(abs(Tq(:)));
c2 = max(aaa_denom)/min(aaa_denom);
if to_save
    nameit = [trial_name,'_the_denoms'];
    saveas(gcf, nameit ,'fig');
    saveas(gcf, nameit,'jpg');
    print('-depsc2',nameit);
end


%% error rates
figure
set(0,'defaultTextInterpreter','latex');
plot(ev_pts, err_opt,'linewidth', 3);
hold on;
plot(ev_pts, err_remez,'--r','linewidth',3.5);
plot(ev_pts, err_aaa,':b','linewidth',3.5);
legend({'Optimization','Remez polynomial','AAA'},'Location','SouthWest','Interpreter','latex')
grid on
set(gca,'FontSize',18)

if to_save
    nameit = [trial_name,'_the_error_rates'];
    saveas(gcf, nameit ,'fig');
    saveas(gcf, nameit,'jpg');
    print('-depsc2',nameit);
end


% printout
 if to_save
     fileID = fopen('summary.txt','w');
     fprintf(fileID,'Sup norm: opt %4.2e AAA %4.2e Remez %4.2e \n', max(abs(err_opt)), max(abs(err_aaa)),  max(abs(err_remez)) );
     fprintf(fileID,'Conditining bound: opt %4.2f AAA %4.2f \n', c1, c2);
     fprintf(fileID,'Mat. Approx. error: opt %4.2e AAA %4.2e Poly Remez %4.2e \n', err_opt_mat, err_aaa_mat, err_poly_mat );
     fclose(fileID);
 else
     fprintf('<strong> Sup norm:</strong> opt %4.2e AAA %4.2e Poly Remez %4.2e \n', max(abs(err_opt)), max(abs(err_aaa)),  max(abs(err_remez)) );
     fprintf('<strong> Conditining bound: </strong> opt %4.2f AAA %4.2f \n', c1, c2);
     fprintf('<strong> Mat. Approx. error:</strong> opt %4.2e AAA %4.2e Poly Remez %4.2e \n', err_opt_mat, err_aaa_mat, err_poly_mat );
 end
    
% save data and close
if to_save    
    save([trial_name,'_data']);
    cd '../'
end



