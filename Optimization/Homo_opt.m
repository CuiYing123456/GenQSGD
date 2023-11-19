clc;clear all;

%% == Start ==
%%%%%%%%%%%% Constant %%%%%%%%%%
N = 10;
alpha = 2e-28;
F_0 = 3e9;
C_0 = 100;  
p_0 = 20;
r_0 = 7.5e7;

D = 128*784;
Q = 2.5;        
G = 33.63;    
L = 0.084;    
sigma = 33.18;   

%% =========== Setting for Step Size Rules: Varying Convergence =============

% ---- System Parameter --------
F_mean = 1e9;
F_ratio = 10;
F_min = 2*F_mean/(1+F_ratio);
F_max = 2*F_ratio*F_mean/(1+F_ratio);
F_n = F_max*ones(floor(N/2), 1);
F_n = [F_n; F_min*ones(ceil(N/2), 1)];
C_n_max = 100e6;
C_n = C_n_max*ones(N, 1);

bits_n = [14];
bits_0 = [32];
s_n = 2^bits_n * ones(N,1);
s_0 = 2^bits_0;

p_mean = 1.5;
r_mean = 50e6 / N;
p_n = p_mean*ones(N,1);
r_n = r_mean*ones(N,1);

% ---------- Setting for Simulation --------
C_max_list = 0.2:0.01:0.25;
T_max = 10000;
gamma_C_list = [0.01];
gamma_E_list = [0.02];
rho_E = 0.9995;
gamma_D_list = [0.02];
rho_D = 600;

%% ========= Simulation 1: Optimize Step Size ========

% --------- Simulation 1 ---------------
energy_line_opt_stepsize = [];
fprintf('Simulation : Optimize Step Size\n');
for i = 1 : length(C_max_list)
    C_max = C_max_list(i);
    fprintf('C_max=%0.3f start ...\n', C_max);
    [ energy_opt_stepsize, K, B, K_0, gamma ] = E_Q_opt_stepsize( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, G, L, sigma, C_max, T_max );
    fprintf('energy_opt=%0.1f\n',  energy_opt_stepsize);
    energy_line_opt_stepsize = [energy_line_opt_stepsize, energy_opt_stepsize];
end

%% ========== Simulation 2: Constant Step Size Rule =============
% ---------- Setting for Simulation --------
gamma_C = 0.01;
% --------- Simulation 2 ---------------
energy_line_cons = [];
K_line_cons = [];    B_line_cons = [];    K_0_line_cons = [];
fprintf(strcat('Simulation : Constant Step Size (gamma=', num2str(gamma_C), ')\n'));
for i = 1 : length(C_max_list)
    C_max = C_max_list(i);
    fprintf('C_max=%0.3f start ...\n', C_max);
    [ energy_cons, K, B, K_0] = E_Q_cons( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_C, G, L, sigma, C_max, T_max );
    fprintf('energy_cons=%0.1f\n',  energy_cons);
    energy_line_cons = [energy_line_cons, energy_cons];
end


%% ========== Simulation 3: Exponential Step Size Rule =============
% ---------- Setting for Simulation --------
gamma_E = 0.02;
% --------- Simulation 3 ---------------
energy_line_exp = [];
fprintf(strcat('Simulation : Exponential Step Size (gamma=', num2str(gamma_E), ')\n'));
for i = 1 : length(C_max_list)
    C_max = C_max_list(i);
    fprintf('C_max=%0.3f start ...\n', C_max);
    [ energy_exp, K, B, K_0 ] = E_Q_exp( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_E, rho_E, G, L, sigma, C_max, T_max );
    fprintf('energy_exp=%0.1f\n',  energy_exp);
    energy_line_exp = [energy_line_exp, energy_exp];
end


%% ========== Simulation 4: Diminishing Step Size Rule =============
% ---------- Setting for Simulation --------
gamma_D = 0.02;
% --------- Simulation 4 ---------------
energy_line_dim = [];
K_line_dim = [];    B_line_dim = [];    K_0_line_dim = [];
fprintf(strcat('Simulation : Diminishing Step Size (gamma=', num2str(gamma_D), ')\n'));
for i = 1 : length(C_max_list)
    C_max = C_max_list(i);
    fprintf('C_max=%0.3f start ...\n', C_max);
    [ energy_dim, K, B, K_0 ] = E_Q_dim( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_D, rho_D, G, L, sigma, C_max, T_max );
    fprintf('energy_dim=%0.1f\n',  energy_dim);
    energy_line_dim = [energy_line_dim, energy_dim];
end

