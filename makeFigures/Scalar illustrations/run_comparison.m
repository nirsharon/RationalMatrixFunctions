function [] = run_comparison(fun, a, b, n, m, l, u, to_save, trial_name)
% A generic function for comparison
%
% rational approx (n,m) over [a,b] for 'fun'
% 
% NS,EK June 2021


% number of coeffs 
n_coefs = n+1;  % numerator
m_coefs = m+1;  % denominator

% sampling
N   = 400; % 512;
pts = linspace(a, b, N);
pts = pts(:);

% N = 300;
%  pts = zeros(N,1);
%  for k=0:N-1 
%      pts(k+1)=.5*(a+b)+.5*(b-a)*cos(pi*(k+0.5)/N); 
%  end


% error evaluation
ev_N = 10^3;
ev_pts = linspace(a, b, ev_N);
ev_pts = ev_pts(:);

% our approximation
tic; 
eps = 1e-15;
[p, q, ~] = RationalMinMaxOpt(fun, n_coefs, m_coefs, pts, l, u, a, b, eps); 
time_opt = toc;

% evaluate the result
p(1) = 2*p(1);
q(1) = 2*q(1);
Tp   = chebeval_scalars(p, ev_pts, n_coefs, a, b);
Tq   = chebeval_scalars(q, ev_pts, m_coefs, a, b);
app  = Tp(:)./Tq(:);
err_opt = (app - fun(ev_pts));

% Remez from "chebfun" package
f1 = chebfun(fun,[a,b]);
tic
[~, q_poly, r] = minimax(f1, n, m); 
time_remez = toc;
err_remez = f1(ev_pts)-r(ev_pts);

% AAA 
y = fun(pts);
tic
max_iters = 100;
[r1, pol, res, zer, z, f, w, errvec]=aaa(y,pts,'degree',n,'lawson',max_iters);
time_aaa = toc;
err_aaa = f1(ev_pts)-r1(ev_pts);

% denominator of AAA
aaa_pts = setdiff(ev_pts, z);
cauchy_mat = 1./bsxfun(@minus,aaa_pts,z.'); 
aaa_denom  = abs(cauchy_mat*w.*(prod(aaa_pts'-z)'));

% open a folder if we need to save
if to_save
    folder_name = [trial_name,'_',datestr(now,'mmmm_dd_yy')];
    mkdir(folder_name)
    cd(folder_name)
end
    
%% the function
figure
set(0,'defaultTextInterpreter','latex');
plot(ev_pts, fun(ev_pts),'linewidth', 3);
set(gca,'FontSize',18)

if to_save
    nameit = [trial_name,'_the_function'];
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
legend('Optimization','Remez','AAA')
grid on
set(gca,'FontSize',18)
e1 = max(abs(err_opt));
e2 = max(abs(err_aaa));
e3 = max(abs(err_remez));
if to_save
    nameit = [trial_name,'_error_rates'];
    saveas(gcf, nameit ,'fig');
    saveas(gcf, nameit,'jpg');
    print('-depsc2',nameit);
end

%% denominaor values
figure
semilogy(ev_pts, Tq(:),'linewidth', 3);
hold on;
semilogy(ev_pts, abs(q_poly(ev_pts)),'--r','linewidth',3.5);
semilogy(aaa_pts, aaa_denom,':b','linewidth',3.5);
legend('Optimization','Remez','AAA','Location','SouthEast')
grid on
set(gca,'FontSize',18)
c1 = max(abs(Tq(:)))/min(abs(Tq(:)));
c2 = max(aaa_denom)/min(aaa_denom);
c3 = max(abs(q_poly(ev_pts)))/min(abs(q_poly(ev_pts)));
if to_save
    nameit = [trial_name,'_the_denoms'];
    saveas(gcf, nameit ,'fig');
    saveas(gcf, nameit,'jpg');
    print('-depsc2',nameit);
end

% printout
if to_save
    fileID = fopen('summary.txt','w');
    fprintf(fileID,'Sup norm: opt %4.2e AAA %4.2e Remez %4.2e \n', e1, e2 ,e3);
    fprintf(fileID,'Timing: opt %4.2f AAA %4.2f Remez %4.2f \n', time_opt, time_aaa, time_remez);
    fprintf(fileID,'Conditining bound: opt %4.2f AAA %4.2f Remez %4.2f \n', c1, c2 ,c3);
    fclose(fileID);
else
    fprintf('<strong> Sup norm:</strong> opt %4.2e AAA %4.2e Remez %4.2e \n', e1, e2 ,e3)
  %  fprintf('Timing: opt %4.2f AAA %4.2f Remez %4.2f \n', time_opt, time_aaa, time_remez)
    fprintf('<strong> Conditining bound: </strong> opt %4.2f AAA %4.2f Remez %4.2f \n', c1, c2 ,c3)
end
    
% save data and close
if to_save    
    save([trial_name,'_data']);
    cd '../'
end


end

