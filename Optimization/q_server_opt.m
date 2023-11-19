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
s_n = 2^32 * ones(N,1);

p_mean = 1.5;
r_mean = 50e6 / N;
p_n = p_mean*ones(N,1);
r_n = r_mean*ones(N,1);

% ---------- Setting for Simulation --------
T_max = 100000;
C_max = 0.25;
bits_list = [10:16, 20, 24, 28, 32];

%% ========= Simulation 1: Proposed Exponential ========
energy_line_Gen = [];
fprintf('Simulation : GenQSGD Opt\n');
for i = 1 : length(bits_list)
    s_0 = 2^bits_list(i);
    fprintf('s_n=%d start ...\n', s_0);
    [ energy_Gen ] = E_Q_opt_stepsize( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, G, L, sigma, C_max, T_max );
    fprintf('energy_Gen=%0.1f\n',  energy_Gen);
    energy_line_Gen = [energy_line_Gen, energy_Gen];
end
