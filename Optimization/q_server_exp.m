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
%% ============ Setting for Step Size Rules: Varying Computation Heterogeneous =============
F_mean = 1e9;
F_ratio = 10;
F_min = 2*F_mean/(1+F_ratio);
F_max = 2*F_ratio*F_mean/(1+F_ratio);
F_n = F_max*ones(floor(N/2), 1);
F_n = [F_n; F_min*ones(ceil(N/2), 1)];
C_n_max = 100e6;    % 100e6
C_n = C_n_max*ones(N, 1);
s_n = 2^32 * ones(N,1);

p_mean = 1.5;
r_mean = 50e6 / N;
p_n = p_mean*ones(N,1);
r_n = r_mean*ones(N,1);

% ---------- Setting for Simulation --------
T_max = 100000;
C_max = 0.25;
bits_list = [10:16, 20, 24, 28, 32];
gamma_exp = 0.02;
rho_E = 0.9995;

%% ========= Simulation 1: Proposed Exponential ========
energy_line_exp = [];
fprintf('Simulation : Proposed Exp\n');
for i = 1 : length(bits_list)
    s_0 = 2^bits_list(i);
    fprintf('s_n=%d start ...\n', s_0);
    [ energy_exp ] = E_Q_exp( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_exp,rho_E, G, L, sigma, C_max, T_max );
    fprintf('energy_opt=%0.1f\n',  energy_exp);
    energy_line_exp = [energy_line_exp, energy_exp];
end

%% ========== Simulation 2: PM-SGD =============
energy_line_PM = [];
fprintf(strcat('Simulation : PM-SGD-opt\n'));
for i = 1 : length(bits_list)
    s_0 = 2^bits_list(i);
    fprintf('s_n=%d start ...\n', s_0);
    [ energy_PM] = E_Q_PM( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_exp, rho_E, G, L, sigma, C_max, T_max );
    fprintf('energy_PM=%0.1f\n',  energy_PM);
    energy_line_PM = [energy_line_PM, energy_PM];
end

%% ========== Simulation 3: FedAvg =============
% --------- Simulation 3 ---------------
energy_line_FedAvg = [];
fprintf(strcat('Simulation : FedAvg-opt\n'));
for i = 1 : length(bits_list)
    s_0 = 2^bits_list(i);
    fprintf('s_n=%d start ...\n', s_0);
    [ energy_FedAvg] = E_Q_FedAvg( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_exp, rho_E, G, L, sigma, C_max, T_max );
    fprintf('energy_FedAvg=%0.1f\n',  energy_FedAvg);
    energy_line_FedAvg = [energy_line_FedAvg, energy_FedAvg];
end

%% ========== Simulation 4: PR-SGD =============
energy_line_PR = [];
fprintf(strcat('Simulation : PR-SGD-opt\n'));
for i = 1 : length(bits_list)
    s_0 = 2^bits_list(i);
    fprintf('s_n=%d start ...\n', s_0);
    [ energy_PR ] = E_Q_PR( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_exp, rho_E, G, L, sigma, C_max, T_max );
    fprintf('energy_PR=%0.1f\n',  energy_PR);
    energy_line_PR = [energy_line_PR, energy_PR];
end


