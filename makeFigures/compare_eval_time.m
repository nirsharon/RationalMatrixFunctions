% script name: "compare_eval_time"
%
% Fixing a matrix A and a vector v, the script compares the timing of
% evaluating F(A)*v via two methods:
% 1 -- rational approximation, r(x) = p(x)/q(x), done by
%      F(A)*v \approx x = r(A)*v, where x solves q(A)x = p(A)v.
% 2 -- Decompose (assuming A is symmetric) A = Q D Q^T. Then, applying
%      (element-wise) f(D) and finally Q*( f(D)* (Q^T*v))) (three mat-vec)

clear; close all;

% matrix size
N_arr = 100:400:1700;

% filter function:
r = .1;    % rise
c = .4;    % center
R = .1;    % Radius
f = @(x) .5*(1 - erf((2/r)*(abs(x-c) - R ) ));
a = -1;    % segment [a,b]
b = 1;

% plot the function:
figure
X = linspace(a,b,200);
plot(X, f(X),'LineWidth',3);

% parameters for the rational app
m = 6;
n = 6;
l = 1;
u = 100;

% sampling points
pts = linspace(a, b, 512);
pts = pts(:);

% calculate the rational approx
[p, q, u] = RationalMinMaxOpt(f, n+1, m+1, pts, l, u, a, b, eps);
p(1) = 2*p(1);
q(1) = 2*q(1);

% plot the approx:
Tp   = chebeval_scalars(p, X, n+1, a, b);
Tq   = chebeval_scalars(q, X, m+1, a, b);
app  = Tp(:)./Tq(:);
hold on
plot(X, app,'LineWidth',3);

e1 = zeros(length(N_arr),1);
e2 = zeros(length(N_arr),1);
t1 = zeros(length(N_arr),1);
t2 = zeros(length(N_arr),1);

for j=1:length(N_arr)
    j
    N = N_arr(j);
    % generate ground truth
    [Q, ~] = qr(randn(N));      % (uniform) random orthogonal
    D = diag(rand(N,1));        % random (uniform) spectrum
    
    % % cluster some parts of the spectrum (optional)
    % k = 5;
    % small_C = 1e-10;
    % D(1:k,1:k) = diag(10*small_C*(1:k)) + D(k+1,k+1);
    
    % the mat and vec
    A = Q*D*Q';
    v = randn(N,1);
    % the groundtruth vector
    FAv = Q*(f(D)*(Q'*v));
    
    
    % METHOD 1: rational approx
    % =====> we assume r is pre-calculated <=====
    tic; % evaluation time for our method
    qA  = chebeval_matrix(q, A, length(q), a, b);
    pAv = cheb_eval_mat_vec(p, A ,m , v);
    FAv_tilde1   = qA\pAv;
    
    % timing
    t1(j) = toc;
    % relative error
    e1(j) = norm(FAv_tilde1 - FAv)/norm(FAv);
    
    %% METHOD 2: spectral decomposition
    tic;
    [V, D_tilde] = eig(A);   % sanity check: norm(V*D_tilde*V' - A)
    FAv_tilde2   = V*(f(D_tilde)*(V'*v));
    
    % timing
    t2(j) = toc;
    % relative error
    e2(j) = norm(FAv_tilde2 - FAv)/norm(FAv);
end

%% summary
figure;
set(0,'defaultTextInterpreter','latex');
plot(N_arr, t1,'LineWidth',3);
hold on
plot(N_arr, t2,'LineWidth',3);
xlabel('Matrix size');
ylabel('Time (seconds)')
legend('Rational','Spectral','Location','NorthWest')
set(gca,'FontSize',18)

