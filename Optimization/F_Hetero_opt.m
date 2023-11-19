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
C_n_max = 100e6;    % 100e6
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
C_max = 0.25;
T_max = 1000000;

%  ------- Setting for System Parameter ---------
F_ratio_list = [ 1:2:10 ];
F_mean = 1e9;

%% ========= Simulation 1: Proposed Exponential ========
energy_line_Gen = [];
fprintf('Simulation : GenQSGD Opt\n');
for i = 1 : length(F_ratio_list)
    F_ratio = F_ratio_list(i);
    F_min = 2*F_mean/(1+F_ratio);
    F_max = 2*F_ratio*F_mean/(1+F_ratio);
    F_n = F_max*ones(floor(N/2), 1);
    F_n = [F_n; F_min*ones(ceil(N/2), 1)];
    fprintf('F_ratio=%d start ...\n', F_ratio);
    [ energy_Gen ] = E_Q_opt_stepsize( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, G, L, sigma, C_max, T_max );
    fprintf('energy_Gen=%0.1f\n',  energy_Gen);
    energy_line_Gen = [energy_line_Gen, energy_Gen];
end

%% ========== Simulation 2: PM-SGD =============
energy_line_PM = [];
fprintf(strcat('Simulation : PM-SGD \n'));
for i = 1 : length(F_ratio_list)
    F_ratio = F_ratio_list(i);
    F_min = 2*F_mean/(1+F_ratio);
    F_max = 2*F_ratio*F_mean/(1+F_ratio);
    F_n = F_max*ones(floor(N/2), 1);
    F_n = [F_n; F_min*ones(ceil(N/2), 1)];
    fprintf('F_ratio=%d start ...\n', F_ratio);
    [ energy_PM] = E_Q_PM( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, G, L, sigma, C_max, T_max );
    fprintf('energy_PM=%0.1f\n',  energy_PM);
    energy_line_PM = [energy_line_PM, energy_PM];
end


%% ========== Simulation 3: FedAvg =============
energy_line_FedAvg = [];
fprintf(strcat('Simulation : FedAvg\n'));
for i = 1 : length(F_ratio_list)
    F_ratio = F_ratio_list(i);
    F_min = 2*F_mean/(1+F_ratio);
    F_max = 2*F_ratio*F_mean/(1+F_ratio);
    F_n = F_max*ones(floor(N/2), 1);
    F_n = [F_n; F_min*ones(ceil(N/2), 1)];
    fprintf('F_ratio=%d start ...\n', F_ratio);
    [ energy_FedAvg] = E_Q_FedAvg( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, G, L, sigma, C_max, T_max );
    fprintf('energy_FedAvg=%0.1f\n',  energy_FedAvg);
    energy_line_FedAvg = [energy_line_FedAvg, energy_FedAvg];
end

%% ========== Simulation 4: PR-SGD =============
energy_line_PR = [];
fprintf(strcat('Simulation : PR-SGD \n'));
for i = 1 : length(F_ratio_list)
    F_ratio = F_ratio_list(i);
    F_min = 2*F_mean/(1+F_ratio);
    F_max = 2*F_ratio*F_mean/(1+F_ratio);
    F_n = F_max*ones(floor(N/2), 1);
    F_n = [F_n; F_min*ones(ceil(N/2), 1)];
    fprintf('F_ratio=%d start ...\n', F_ratio);
    [ energy_PR ] = E_Q_PR( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, G, L, sigma, C_max, T_max );
    fprintf('energy_PR=%0.1f\n',  energy_PR);
    energy_line_PR = [energy_line_PR, energy_PR];
end
