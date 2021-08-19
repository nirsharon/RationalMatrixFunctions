% script name: "figure_filter_vec_timing"
%
% Fixing a matrix A and a vector v, the script compares the timing of
% evaluating F(A)*v via two methods:
% 1 -- rational approximation, r(x) = p(x)/q(x), done by
%      F(A)*v \approx x = r(A)*v, where x solves q(A)x = p(A)v.
% 2 -- Decompose (assuming A is symmetric) A = Q D Q^T. Then, applying
%      (element-wise) f(D) and finally Q*( f(D)* (Q^T*v))) (three mat-vec)
%
% NS, August 21

clear;
close all;
to_save = 1;

%% settings: matrix size and repeating
N_arr = 100:200:2500;
repeat_trial = 10;

% filter function:
r = .1;    % rise
c = .4;    % center
R = .1;    % Radius
f = @(x) .5*(1 - erf((2/r)*(abs(x-c) - R ) ));

% parameters for the rational app
m1 = 5;   n1 = 5;
m2 = 10;  n2 = 10;
l = 1;    u = 1000;
a = -1;   b = 1;

%% sampling points
N   = 500;
pts = linspace(a, b, N);
pts = pts(:);

% calculate the rational approxs
[p1, q1] = RationalMinMaxOpt(f, n1+1, m1+1, pts, l, u, a, b, eps);
p1(1) = 2*p1(1);
q1(1) = 2*q1(1);
[p2, q2] = RationalMinMaxOpt(f, n2+1, m2+1, pts, l, u, a, b, eps);
p2(1) = 2*p2(1);
q2(1) = 2*q2(1);

% plot the approxs
X = linspace(a,b,200);
X = X(:);
Tp1   = chebeval_scalars(p1, X, n1+1, a, b);
Tq1   = chebeval_scalars(q1, X, m1+1, a, b);
app1  = Tp1(:)./Tq1(:);
e1 = max(abs(app1-f(X)));
Tp2   = chebeval_scalars(p2, X, n2+1, a, b);
Tq2   = chebeval_scalars(q2, X, m2+1, a, b);
app2  = Tp2(:)./Tq2(:);
e2 = max(abs(app2-f(X)));


%% plot the function:
figure
set(0,'defaultTextInterpreter','latex');
plot(X, f(X), 'r','LineWidth',3);
set(gca,'FontSize',18)
hold on
plot(X, app1,'b','LineWidth',3);
plot(X, app2,'--b','LineWidth',3);
legend({'The filter function','Rational $(5,5)$',...
    'Rational $(10,10)$'},'Location','NorthWest','Interpreter','latex')
fprintf('<strong> Sup norms:</strong> \n type (5,5) is %4.2e \n type (10,10) is %4.2e \n', e1, e2);

if to_save
    trial_name  = 'filter_mat_vecV2';
    folder_name = [trial_name,'_',datestr(now,'mmmm_dd_yy')];
    % make new folder
    mkdir(folder_name)
    cd(folder_name)
    % save figure
    nameit = [trial_name,'_function'];
    saveas(gcf, nameit ,'fig');
    saveas(gcf, nameit,'jpg');
    print('-depsc2',nameit);
end

%% main loop
% initial variables
err_rat1 = zeros(length(N_arr), repeat_trial);
err_rat2 = zeros(length(N_arr), repeat_trial);
err_spec = zeros(length(N_arr), repeat_trial);
time_rat1 = zeros(length(N_arr), repeat_trial);
time_rat2 = zeros(length(N_arr), repeat_trial);
time_spec = zeros(length(N_arr), repeat_trial);

for trialnum = 1:repeat_trial
    fprintf('We start trial number %d \n',trialnum);
    fprintf('N size');
    for j=1:length(N_arr)
        N = N_arr(j);
        fprintf(',%d ',N);
        % generate ground truth
        [Q, ~] = qr(randn(N));      % (uniform) random orthogonal
        D = diag(-1+2*rand(N,1));   % random (uniform) spectrum
        
        % the mat and vec
        A = Q*D*Q';
        v = randn(N,1);
        
        % the groundtruth vector
        FAv = Q*(f(D)*(Q'*v));
        
        % METHOD 1: rational approx (5,5)
        % =====> we assume r is pre-calculated <=====
        tic; % evaluation time for our method
        qA1  = chebeval_matrix(q1, A, length(q1), a, b);
        pAv1 = cheb_eval_mat_vec(p1, A ,m1 + 1 , a, b, v);
        FAv_tilde1   = qA1\pAv1;
        time_rat1(j,trialnum) = toc;
        err_rat1(j,trialnum)  = norm(FAv_tilde1 - FAv)/norm(FAv);
        
%         % sanity check:
%         pA1=chebeval_matrix(p1,A,length(p1),a,b);
%         norm(FAv_tilde1 - FAv)/norm(v) < norm(Q*f(D)*Q'-qA1\pA1)
        
        % METHOD 2: rational approx (10,10)
        % =====> we assume r is pre-calculated <=====
        tic; % evaluation time for our method
        qA2  = chebeval_matrix(q2, A, length(q2), a, b);
        pAv2 = cheb_eval_mat_vec(p2, A ,m2 + 1 , a, b, v);
        FAv_tilde2   = qA2\pAv2;
        time_rat2(j,trialnum) = toc;
        err_rat2(j,trialnum) = norm(FAv_tilde2 - FAv)/norm(FAv);
        
        %% METHOD 3: spectral decomposition
        tic;
        [V, D_tilde] = eig(A);   % sanity check: norm(V*D_tilde*V' - A)
        FAv_tilde3   = V*(f(D_tilde)*(V'*v));
        time_spec(j,trialnum) = toc;
        err_spec(j,trialnum) = norm(FAv_tilde3 - FAv)/norm(FAv);
    end
    fprintf('. \n')
end

%% summary
figure;
set(0,'defaultTextInterpreter','latex');
plot(N_arr, mean(time_rat1,2),'b','LineWidth',3);
hold on
plot(N_arr, mean(time_rat2,2),'--b','LineWidth',3);
plot(N_arr, mean(time_spec,2),':r','LineWidth',3);
xlabel('Matrix size');
ylabel('Time (seconds)')
legend({'Rational (5,5)','Rational (10,10)','Spectral'},'Location','NorthWest','Interpreter','latex')
set(gca,'FontSize',18)

if to_save
    % save figure
    nameit = [trial_name,'_timing'];
    saveas(gcf, nameit ,'fig');
    saveas(gcf, nameit,'jpg');
    print('-depsc2',nameit);
    % save data
    save([trial_name,'_data']);
    cd '../'
end

figure;
set(0,'defaultTextInterpreter','latex');
semilogy(N_arr, mean(err_rat1,2),'LineWidth',3);
hold on
grid on
semilogy(N_arr, mean(err_rat2,2),'LineWidth',3);
semilogy(N_arr, mean(err_spec,2),'LineWidth',3);
xlabel('Matrix size');
ylabel('Error (seconds)')
legend({'Rational (5,5)','Rational (5,5)','Spectral'},'Location','NorthWest','Interpreter','latex')
set(gca,'FontSize',18)
