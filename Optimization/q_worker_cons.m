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
C_n_max = 100e6;   
C_n = C_n_max*ones(N, 1);
s_0 = 2^32;

p_mean = 1.5;
r_mean = 50e6 / N;
p_n = p_mean*ones(N,1);
r_n = r_mean*ones(N,1);

% ---------- Setting for Simulation --------
T_max = 100000;
C_max = 0.25;
bits_list = [10:16, 20, 24, 28, 32];
gamma_cons = 0.01; 

%% ========= Simulation 1: Proposed Constant ========
energy_line_cons = [];
fprintf('Simulation : Proposed Constant\n');
for i = 1 : length(bits_list)
    s_n = 2^bits_list(i) * ones(N,1);   
    fprintf('s_n=%d start ...\n', s_n(1));
    [ energy_cons ] = E_Q_cons( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_cons, G, L, sigma, C_max, T_max );
    fprintf('energy_opt=%0.1f\n',  energy_cons);
    energy_line_cons = [energy_line_cons, energy_cons];
end

%% ========== Simulation 2-1: PM-SGD =============
energy_line_PM = [];
fprintf(strcat('Simulation : PM-SGD-opt\n'));
for i = 1 : length(bits_list)
    s_n = 2^bits_list(i) * ones(N,1);
    fprintf('s_n=%d start ...\n', s_n(1));
    [ energy_cons] = E_Q_PM( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_cons, G, L, sigma, C_max, T_max );
    fprintf('energy_cons=%0.1f\n',  energy_cons);
    energy_line_PM = [energy_line_PM, energy_cons];
end

%% ========== Simulation 3-1: FedAvg =============
energy_line_FedAvg = [];
fprintf(strcat('Simulation : FedAvg-opt\n'));
for i = 1 : length(bits_list)
    s_n = 2^bits_list(i) * ones(N,1);
    fprintf('s_n=%d start ...\n', s_n(1));
    [ energy_FedAvg] = E_Q_FedAvg( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_cons, G, L, sigma, C_max, T_max );
    fprintf('energy_FedAvg=%0.1f\n',  energy_FedAvg);
    energy_line_FedAvg = [energy_line_FedAvg, energy_FedAvg];
end

%% ========== Simulation 4-1: PR-SGD =============
energy_line_PR = [];
fprintf(strcat('Simulation : PR-SGD-opt\n'));
for i = 1 : length(bits_list)
    s_n = 2^bits_list(i) * ones(N,1);
    fprintf('s_n=%d start ...\n', s_n(1));
    [ energy_PR ] = E_Q_PR( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_cons, G, L, sigma, C_max, T_max );
    fprintf('energy_PR=%0.1f\n',  energy_PR);
    energy_line_PR = [energy_line_PR, energy_PR];
end

