clc; clear; close all;

%% Parameters
L  = 1;
T  = 1;

h = 0.02;    %L / Nx;
tau = 0.0002; %T / Nt;
r  = tau / h^2;
Nx = round(L/h);                 % number of spatial intervals
Nt = round(T/tau);               % number of time steps

x = linspace(0,L,Nx+1)';
t = linspace(0,T,Nt+1);
Q = @(x,t) (4*pi^2 - 1) * exp(-t) .* sin(2*pi*x); % Source term



%%  Part1: solution at final time

% Exact solution
u_exact = @(x,t) exp(-t) .* sin(2*pi*x);

%%Initial condition
u0 = sin(2*pi*x);
% ===============================
u_ex_T = u_exact(x, T);
U_exp = euler_explicit(Nx,Nt,u0, r, tau, x, t,Q);
U_cn = crank_nicolson(Nx,Nt,u0, r, tau, x, t,Q);
U_imp = euler_implicit(Nx,Nt,u0, r, tau, x, t,Q);
% ===============================
%  Plots: Solution at final time
%===============================
figure;
plot(x, u_ex_T, 'k-', 'LineWidth', 1); hold on;
plot(x, U_exp(:,end), 'r--','LineWidth', 4);  hold on;
plot(x, U_imp(:,end), 'b-.','LineWidth', 2);
plot(x, U_cn(:,end),  'g:','LineWidth', 2);
legend('Exact','Explicit','Implicit','Crank-Nicolson','Location','best');
xlabel('x'); ylabel('u(x,T)');
title('Solution at final time T=1; h=0.02 and tau=0.0002');
grid on;

% Error plots
% ===============================
err_exp = abs(U_exp(:,end) - u_ex_T);
err_imp = abs(U_imp(:,end) - u_ex_T);
err_cn  = abs(U_cn(:,end)  - u_ex_T);

figure;
plot(x, err_exp, 'r--', 'LineWidth', 1.5); hold on;
plot(x, err_imp, 'b-.', 'LineWidth', 1.5);
plot(x, err_cn,  'g:', 'LineWidth', 1.5);
legend('Explicit','Implicit','Crank-Nicolson','Location','best');
xlabel('x'); ylabel('Absolute Error');
title('Error at final time T=1');
grid on;

figure;
mesh(t, x, U_exp) % or U_imp, U_cn
xlabel('Time t'); ylabel('x'); zlabel('u(x,t)');
title('Evolution of u(x,t) over time');


%% Part 2: Eigen Values and Condition Number
h_values = [0.2, 0.1, 0.05, 0.025];
K_values = zeros(size(h_values));
lambda_max = zeros(size(h_values));
lambda_min = zeros(size(h_values));
for k = 1:length(h_values)
    dx = h_values(k);
    Nx = L/dx;
    A_spatial = spdiags([ones(Nx-1,1), -2*ones(Nx-1,1), ones(Nx-1,1)], -1:1, Nx-1, Nx-1) / dx^2;
    K_values(k) = cond(A_spatial,2);
    lambda = eig(A_spatial);
    lambda_max(k) = max(lambda);
    lambda_min(k) = min(lambda);
    fprintf('Nx=%d, dx=%.4f, lambda_max=%.4f, lambda_min=%.4f, K2=%.4f\n', Nx, dx, lambda_max(k), lambda_min(k), K_values(k));

end

figure;
plot(h_values, K_values, 'o-');
xlabel('h (dx)'); ylabel('Condition number K2(A)');
title('Dependence of K2(A) on spatial step-size h');
grid on;

%% Part 3: Stability

% (h, tau) cases to test
h_vals   = [0.1  0.1 0.02];
tau_vals = [0.01 0.0025 0.0002];
L  = 1;
T  = 1;
Uexp_store = cell(length(h_vals),1);
x_store    = cell(length(h_vals),1);

for k = 1:length(h_vals)

    h   = h_vals(k);
    tau = tau_vals(k);

    Nx_2 = L / h;
    Nt_2 = T / tau;
    r_2  = tau / h^2;

    x_2 = linspace(0,L,Nx_2+1)';
    t_2 = linspace(0,T,Nt_2+1);

    u0_2 = sin(2*pi*x_2);

    fprintf('\nCase %d: h = %.3f, tau = %.4f, r = %.3f\n', k, h, tau, r_2);

    if r > 0.5
        fprintf('  -> Explicit Euler expected to be UNSTABLE\n');
    else
        fprintf('  -> Explicit Euler expected to be STABLE\n');
    end    
    
    U_exp_2 = euler_explicit(Nx_2, Nt_2, u0_2, r_2, tau, x_2, t_2, Q);
    U_imp_2 = euler_implicit(Nx_2, Nt_2, u0_2, r_2, tau, x_2, t_2, Q);
    U_cn_2  = crank_nicolson(Nx_2, Nt_2, u0_2, r_2, tau, x_2, t_2, Q);
    Uexp_store{k} = U_exp_2(:,end);
    x_store{k}    = x;

    figure
    sgtitle(sprintf('h = %.2f, tau = %.4f', h, tau))
    plot(x_2, U_exp_2(:,end),'r--', 'LineWidth', 1.5)
    hold on
    plot(x_2, U_imp_2(:,end),'b-.', 'LineWidth', 1.5)
    plot(x_2, U_cn_2(:,end), 'g:','LineWidth', 1.5)
    hold off

    legend('Explicit','Implicit','Crank-Nicolson','Location','best')
    xlabel('x'); ylabel('u(x,T)')
    grid on

end
% figure
% hold on
% for k = 1:length(h_vals)
%     plot(x_store{k}, Uexp_store{k}, 'LineWidth', 1.5, ...
%         'DisplayName', sprintf('h=%.3f, \\tau=%.4f', h_vals(k), tau_vals(k)));
% end
% hold off
% grid on
% xlabel('x')
% ylabel('u(x,T)')
% title('Explicit Euler solution at T=1 for different (h, \tau)')
% legend('Location','best')


%%
L = 1; T = 1;
h0 = 0.1;    %L / Nx;
tau0 = 0.01; %T / Nt;

%u0 = sin(2*pi*x);
levels = 3;

err_exp = zeros(levels,1);
err_imp = zeros(levels,1);
err_cn  = zeros(levels,1);

for k = 1:levels

    h   = h0   / 2^(k-1);
    tau = tau0 / 2^(k-1);
    fprintf('h = %.6f\n', h)
    fprintf('tau = %.4f\n', tau)
    Nx_3 = L / h;
    Nt_3 = T / tau;
    r  = tau / h^2;

    x_3 = linspace(0,L,Nx_3+1)';
    t_3 = linspace(0,T,Nt_3+1);

    u0_3 = sin(2*pi*x_3);

    % Numerical solutions
    U_exp_3 = euler_explicit(Nx_3,Nt_3,u0_3,r,tau,x_3,t_3,Q);
    U_imp_3 = euler_implicit(Nx_3,Nt_3,u0_3,r,tau,x_3,t_3,Q);
    U_cn_3  = crank_nicolson(Nx_3,Nt_3,u0_3,r,tau,x_3,t_3,Q);

    % Exact solution on grid
    [X,Tm] = meshgrid(x_3,t_3);
    U_ex_3 = exp(-Tm').*sin(2*pi*X');
    

    % L∞(0,1;L∞(0,1)) error
    err_exp(k) = max(max(abs(U_exp_3 - U_ex_3)));
    err_imp(k) = max(max(abs(U_imp_3 - U_ex_3)));
    err_cn(k)  = max(max(abs(U_cn_3  - U_ex_3)));

end

% Observed convergence orders
p_exp = log(err_exp(1:end-1)./err_exp(2:end)) / log(2);
p_imp = log(err_imp(1:end-1)./err_imp(2:end)) / log(2);
p_cn  = log(err_cn (1:end-1)./err_cn (2:end)) / log(2);

disp('Observed orders:')
disp(table(p_exp, p_imp, p_cn, ...
    'VariableNames',{'Explicit','Implicit','CrankNicolson'}))

hvals = [h0, h0/2, h0/4];

figure
loglog(hvals, err_imp, 'o-', 'LineWidth', 1.5); hold on
loglog(hvals, err_cn,  's-', 'LineWidth', 1.5)
grid on
xlabel('h')
ylabel('L^\infty error')
legend('Implicit Euler','Crank–Nicolson','Location','best')

%%
figure;
mesh(t_3, x_3, U_ex_3) % or U_imp, U_cn
xlabel('Time t'); ylabel('x'); zlabel('u(x,t)');
title('Evolution of u(x,t) over time');

%% Functions
%% ===============================
%  Explicit Euler (FTCS)
function [U_explicit] = euler_explicit(Nx,Nt,u0, r, tau, x, t,Q)
U_explicit = zeros(Nx+1, Nt+1);
U_explicit(:,1) = u0;
for n = 1:Nt
    for i = 2:Nx
        U_explicit(i,n+1) = U_explicit(i,n) ...
            + r*(U_explicit(i+1,n) - 2*U_explicit(i,n) + U_explicit(i-1,n)) ...
            + tau * Q(x(i), t(n));
    end
end
end

%% ===============================
%  Implicit Euler (BTCS)
function [U_implicit] = euler_implicit(Nx,Nt,u0, r, tau, x, t,Q)
U_implicit = zeros(Nx+1, Nt+1);
U_implicit(:,1) = u0;
A = spdiags([-r*ones(Nx-1,1), (1+2*r)*ones(Nx-1,1), -r*ones(Nx-1,1)], -1:1, Nx-1, Nx-1);

for n = 1:Nt
    b = U_implicit(2:Nx,n) + tau * Q(x(2:Nx), t(n+1));
    U_implicit(2:Nx,n+1) = A \ b;
end
end
%% ===============================
%  Crank–Nicolson
function [U_crank] = crank_nicolson(Nx,Nt,u0, r, tau, x, t,Q)
U_crank  = zeros(Nx+1, Nt+1);
U_crank(:,1)  = u0;
A_cn = spdiags([-r/2*ones(Nx-1,1), (1+r)*ones(Nx-1,1), -r/2*ones(Nx-1,1)], -1:1, Nx-1, Nx-1);
B_cn = spdiags([ r/2*ones(Nx-1,1), (1-r)*ones(Nx-1,1),  r/2*ones(Nx-1,1)], -1:1, Nx-1, Nx-1);

for n = 1:Nt
    q = tau/2 * ( Q(x(2:Nx), t(n)) + Q(x(2:Nx), t(n+1)) );
    rhs = B_cn * U_crank(2:Nx,n) + q;
    U_crank(2:Nx,n+1) = A_cn \ rhs;
end
end
