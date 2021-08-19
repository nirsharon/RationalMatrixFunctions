% script name: "figure_spd_timing"
%
% Fixing a matrix A and a vector v, the script compares the timing of
% evaluating P_spd(A)*v via two methods:
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
minN = 100;
maxN = 2500;
Nlen = 13;
N_arr = floor(linspace(minN, maxN, Nlen));
repeat_trial = 10;

if to_save
    trial_name  = 'SPD_proj_vec';
    folder_name = [trial_name,'_',datestr(now,'mmmm_dd_yy')];
    % make new folder
    mkdir(folder_name)
    cd(folder_name)
end

%% projection function
f = @(T) max(0,T);

% parameters for the rational app
m = 5;   n = 5;

l = 1;    u = 100;
a = -1;   b = 1;

%% sampling points
N   = 500;
pts = linspace(a, b, N);
pts = pts(:);

% calculate the rational approxs
[p1, q1] = RationalMinMaxOpt_pos(f, n+1, m+1, pts, l, u, a, b, eps);
p1(1) = 2*p1(1);
q1(1) = 2*q1(1);

% plot the approxs
X = linspace(a,b,200);
X = X(:);
Tp1   = chebeval_scalars(p1, X, n+1, a, b);
Tq1   = chebeval_scalars(q1, X, m+1, a, b);
app1  = Tp1(:)./Tq1(:);
e1 = max(abs(app1-f(X)));


%% main loop
% initial variables
err_rat1 = zeros(length(N_arr), repeat_trial);
err_rat2 = zeros(length(N_arr), repeat_trial);
err_spec = zeros(length(N_arr), repeat_trial);
err_spec2 = zeros(length(N_arr), repeat_trial);

time_rat1 = zeros(length(N_arr), repeat_trial);
time_rat2 = zeros(length(N_arr), repeat_trial);
time_spec = zeros(length(N_arr), repeat_trial);
time_spec2 = zeros(length(N_arr), repeat_trial);

for trialnum = 1:repeat_trial
    fprintf('We start trial number %d \n',trialnum);
    fprintf('N size');
    for j=1:length(N_arr)
        N = N_arr(j);
        fprintf(',%d ',N);
        
        % generate same orthogonal matrix
        [Q, ~] = qr(randn(N));      % (uniform) random orthogonal
        
        %==========> uniform spectrum

        D1 = diag(-1+2*rand(N,1));   % random (uniform) spectrum
        
        % the mat and vec
        A1 = Q*D1*Q';
        v = randn(N,1);
        
        % the groundtruth vector
        FA1v = Q*(f(D1)*(Q'*v));
        
        % METHOD 1: rational approx (5,5)
        tic;
        qA1  = chebeval_matrix(q1, A1, length(q1), a, b);
        pA1v = cheb_eval_mat_vec(p1, A1 ,m + 1 , a, b, v);
        FAv_tilde1   = qA1\pA1v;
        time_rat1(j,trialnum) = toc;
        err_rat1(j,trialnum)  = norm(FAv_tilde1 - FA1v)/norm(FA1v);
        
%         % sanity check:
%         pA1=chebeval_matrix(p1,A,length(p1),a,b);
%         norm(FAv_tilde1 - FA1v)/norm(v) < norm(Q*f(D1)*Q'-qA1\pA1)
        
        % METHOD 2: spectral decomposition
        tic;
        [V, D_tilde] = eig(A1);   % sanity check: norm(V*D_tilde*V' - A)
        FAv_spec   = V*(f(D_tilde)*(V'*v));
        time_spec(j,trialnum) = toc;
        err_spec(j,trialnum) = norm(FAv_spec - FA1v)/norm(FA1v);
        
        %==========> clustered spectrum
        
        % cluster some parts of the spectrum (optional)
        D2 = diag(-1+2*rand(N,1));   % random (uniform) spectrum
        k = floor(N/2);
        small_C = 1e-10;
        D2(1:k,1:k) = diag(0.3 + small_C*(1:k));
        
        D2(k+1:end,k+1:end) = diag(-0.3-small_C*(1:(N-k)));
                
        % the mat and vec
        A2 = Q*D2*Q';
               
        % the groundtruth vector
        FA2v = Q*(f(D2)*(Q'*v));
        
        % METHOD 1: rational approx (5,5)
        tic;
        qA2  = chebeval_matrix(q1, A2, length(q1), a, b);
        pA2v = cheb_eval_mat_vec(p1, A2 ,m + 1 , a, b, v);
        FAv_tilde2   = qA2\pA2v;
        time_rat2(j,trialnum) = toc;
        err_rat2(j,trialnum)  = norm(FAv_tilde2 - FA2v)/norm(FA2v);
                
        % METHOD 2: spectral decomposition
        tic;
        [V, D_tilde] = eig(A2);   % sanity check: norm(V*D_tilde*V' - A)
        FAv_spec2   = V*(f(D_tilde)*(V'*v));
        time_spec2(j,trialnum) = toc;
        err_spec2(j,trialnum) = norm(FAv_spec2 - FA2v)/norm(FA2v);
        
      
    end
    fprintf('. \n')
end

%% summary
fprintf('<strong> Difference in runtime between two cases :</strong> %2.2e \n', norm(mean(time_rat1,2) - mean(time_rat2,2)))

figure;
set(0,'defaultTextInterpreter','latex');
plot(N_arr, mean(time_rat1,2),'b','LineWidth',3);
hold on
%plot(N_arr, mean(time_rat2,2),'--b','LineWidth',3);
plot(N_arr, mean(time_spec,2),':r','LineWidth',3);
plot(N_arr, mean(time_spec2,2),'r','LineWidth',3);
xlabel('Matrix size');
ylabel('Time (seconds)')
legend({'Rational (both cases)', 'Spectral (uniform spectrum)','Spectral (clustered spectrum)' ...
    },'Location','NorthWest','Interpreter','latex')
%legend({'Rational 1','Rational 2', 'Spectral','Spectral2'},'Location','NorthWest','Interpreter','latex')
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


% figure;
% set(0,'defaultTextInterpreter','latex');
% semilogy(N_arr, mean(err_rat1,2),'LineWidth',3);
% hold on
% grid on
% semilogy(N_arr, mean(err_rat2,2),'LineWidth',3);
% semilogy(N_arr, mean(err_spec,2),'LineWidth',3);
% semilogy(N_arr, mean(err_spec2,2),'LineWidth',3);
% xlabel('Matrix size');
% ylabel('Error')
% legend({'Rational 1','Rational 2', 'Spectral','Spectral2'},'Location','NorthWest','Interpreter','latex')
% set(gca,'FontSize',18)

% % plot the histogram
% hist_bins = 100;
% bin_edges = linspace(-1,1,hist_bins); 
% set(0,'defaultTextInterpreter','latex');
% figure; h1 = histogram(diag(D1),bin_edges); 
% hold on; h2 = histogram(diag(D2),bin_edges);
% legend({'D_1','D_2'})
% set(gca,'FontSize',18)
% 
% if to_save
%     nameit = [trial_name,'_hists'];
%     saveas(gcf, nameit ,'fig');
%     saveas(gcf, nameit,'jpg');
%     print('-depsc2',nameit);
% end



