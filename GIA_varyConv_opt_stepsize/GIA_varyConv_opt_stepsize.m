clc;clear all;

%% == Start ==
%%%%%%%%%%%% Constant %%%%%%%%%%
N = 10;
alpha = 2e-28;
beta = 1e-28;
F_0 = 3e9;
C_0 = 100;       % 1000
p_0 = 20;
r_0 = 7.5e7;

D = 128*784;
Q = 2.5;          % 3.2-0.7=2.5
G = 33.63;      % 0.6
L = 0.084;     %1.4e3~1.5e3        %2.4163     22   L=1
sigma = 33.18;     % 18

clr=[[1 0 0];[1 0 1];[0 0 1];[0 0.5 0.5];[0.2 0.5 0.3];[0 0 0]];
str=[' -';' :';'--'; '-.'];
mkr = ['o';'+';'s';'^';'*';'d'];
pathname_tmp = './Tmp/';
pathname_data = './Data/';
rand_stop = 20;

%% =========== Setting for Step Size Rules: Varying Convergence =============
%  ------- Setting for System Parameter ---------
% F_ratio = 10;
% C_ratio = 1;

% ---- System Parameter --------
% F_mean = 2e9;
% F_ratio = 1;
F_mean = 1e9;
F_ratio = 10;
F_min = 2*F_mean/(1+F_ratio);
F_max = 2*F_ratio*F_mean/(1+F_ratio);
F_n = F_max*ones(floor(N/2), 1);
F_n = [F_n; F_min*ones(ceil(N/2), 1)];
% F_n = 2e9*ones(N, 1);
C_n_max = 100e6;    % 100e6
C_n = C_n_max*ones(N, 1);

bits_n = [14];
bits_0 = [32];
s_n = 2^bits_n * ones(N,1);
s_0 = 2^bits_0;
% bits_n = [14];
% bits_0 = [32];
% sn_mean = 2^bits_n;
% sn_ratio = 10;
% sn_min = 2*sn_mean/(1+sn_ratio);
% sn_max = 2*sn_ratio*sn_mean/(1+sn_ratio);
% s_n = sn_max*ones(floor(N/2), 1);
% s_n = [s_n; sn_min*ones(ceil(N/2), 1)];
% s_0 = 2^bits_0;

p_mean = 1.5;
r_mean = 50e6 / N;
p_n = p_mean*ones(N,1);
r_n = r_mean*ones(N,1);

% ---------- Setting for Simulation --------
% C_max_list = 2.8:0.1:3.5;
C_max_list = 0.2:0.01:0.25;
T_max = 10000;
gamma_C_list = [0.01];
gamma_E_list = [0.02];
rho_E = 0.9995;
gamma_D_list = [0.02];
rho_D = 600;

% ----------- Setting for Plot ---------------
figure(1);
Lgd_cnt =1;
energy_max = 0;
energy_min = 1000000;

%% ========= Simulation 1: Optimize Step Size ========

% --------- Simulation 1 ---------------
energy_line_opt_stepsize = [];
K_line_opt = [];    B_line_opt = [];    K_0_line_opt = [];  gamma_line_opt = [];
fprintf('Simulation : Optimize Step Size\n');
for i = 1 : length(C_max_list)
    C_max = C_max_list(i);
    fprintf('C_max=%0.3f start ...\n', C_max);
    [ energy_opt_stepsize, K, B, K_0, gamma ] = E_Q_opt_stepsize( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, G, L, sigma, C_max, T_max );
    fprintf('energy_opt=%0.1f\n',  energy_opt_stepsize);
    energy_line_opt_stepsize = [energy_line_opt_stepsize, energy_opt_stepsize];
    K_line_opt = [K_line_opt, K];   B_line_opt = [B_line_opt, B];   K_0_line_opt = [K_0_line_opt, K_0]; gamma_line_opt = [gamma_line_opt, gamma];
end
save(strcat(pathname_data, 'Energy_opt_stepsize.mat'), 'energy_line_opt_stepsize');
% save(strcat(pathname_data, 'K_line.mat'), 'K_line');
% save(strcat(pathname_data, 'B_line.mat'), 'B_line');
% save(strcat(pathname_data, 'K_0_line.mat'), 'K_0_line');
% save(strcat(pathname_data, 'gamma_line.mat'), 'gamma_line');

% ------------ Plot 1 ---------
p1 = plot(C_max_list, energy_line_opt_stepsize);
%         p2 = loglog(C_max_list, energy_line_opt_stepsize);
p1.Color = clr(1, :);
p1.LineStyle = str(1, :);
p1.Marker = mkr(1, :);
p1.MarkerSize = 8;
p1.LineWidth = 1.5;
Lgd{Lgd_cnt} = 'Opt';
Lgd_cnt = Lgd_cnt + 1;
hold on

energy_max = max([ energy_line_opt_stepsize(energy_line_opt_stepsize~=Inf) ]);
energy_min = min([energy_line_opt_stepsize]);


%% ========== Simulation 2: Constant Step Size Rule =============
% ---------- Setting for Simulation --------
% gamma_C_list = [0.01];
for j = 1 : length(gamma_C_list)
    gamma_C = gamma_C_list(j);
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
        K_line_cons = [K_line_cons, K];   B_line_cons = [B_line_cons, B];   K_0_line_cons = [K_0_line_cons, K_0];
    end
    save(strcat(pathname_data, 'Energy_cons_gamma=', num2str(gamma_C_list(j)), '.mat'), 'energy_line_cons');
    
    % ------------ Plot i ---------
    p2 = plot(C_max_list, energy_line_cons);
    %         p1 = loglog(C_max_list, energy_line_cons);
    p2.Color = clr(2, :);
    p2.LineStyle = str(2, :);
    p2.Marker = mkr(j+1, :);
    p2.MarkerSize = 8;
    p2.LineWidth = 1.5;
	Lgd{Lgd_cnt} = strcat( 'C gamma=', num2str(gamma_C_list(j)) );
    Lgd_cnt = Lgd_cnt + 1;
    hold on
    
    energy_max = max([ energy_line_cons(energy_line_cons~=Inf), energy_max ]);
    energy_min = min([ energy_line_cons, energy_min ]);

end

%% ========== Simulation 3: Exponential Step Size Rule =============
% ---------- Setting for Simulation --------
% gamma_E_list = [0.02];
% rho_E = 0.9995;
for j = 1 : length(gamma_E_list)
    gamma_E = gamma_E_list(j);
    % --------- Simulation 3 ---------------
    energy_line_exp = [];
    K_line_exp = [];    B_line_exp = [];    K_0_line_exp = [];
    fprintf(strcat('Simulation : Exponential Step Size (gamma=', num2str(gamma_E), ')\n'));
    for i = 1 : length(C_max_list)
        C_max = C_max_list(i);
        fprintf('C_max=%0.3f start ...\n', C_max);
        [ energy_exp, K, B, K_0 ] = E_Q_exp( N, s_n, s_0, D, alpha, ...
            C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
            Q, gamma_E, rho_E, G, L, sigma, C_max, T_max );
        fprintf('energy_exp=%0.1f\n',  energy_exp);
        energy_line_exp = [energy_line_exp, energy_exp];
        K_line_exp = [K_line_exp, K];   B_line_exp = [B_line_exp, B];   K_0_line_exp = [K_0_line_exp, K_0];
    end
    save(strcat(pathname_data, 'Energy_exp_gamma=', num2str(gamma_E_list(j)), '.mat'), 'energy_line_exp');
    
    % ------------ Plot i ---------
    p3 = plot(C_max_list, energy_line_exp);
    %         p1 = loglog(C_max_list, energy_line_cons);
    p3.Color = clr(3, :);
    p3.LineStyle = str(3, :);
    p3.Marker = mkr(j+1, :);
    p3.MarkerSize = 8;
    p3.LineWidth = 1.5;
	Lgd{Lgd_cnt} = strcat( 'E gamma=', num2str(gamma_E_list(j)) );
    Lgd_cnt = Lgd_cnt + 1;
    hold on
    
    energy_max = max([ energy_line_exp(energy_line_exp~=Inf), energy_max ]);
    energy_min = min([ energy_line_exp, energy_min ]);

end

%% ========== Simulation 4: Diminishing Step Size Rule =============
% ---------- Setting for Simulation --------
% gamma_D_list = [0.02];
% rho_D = 600;
for j = 1 : length(gamma_D_list)
    gamma_D = gamma_D_list(j);
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
        K_line_dim = [K_line_dim, K];   B_line_dim = [B_line_dim, B];   K_0_line_dim = [K_0_line_dim, K_0];
    end
    save(strcat(pathname_data, 'Energy_dim_gamma=', num2str(gamma_D_list(j)), '.mat'), 'energy_line_dim');
    
    % ------------ Plot i ---------
    p4 = plot(C_max_list, energy_line_dim);
    %         p4 = loglog(C_max_list, energy_line_cons);
    p4.Color = clr(4, :);
    p4.LineStyle = str(4, :);
    p4.Marker = mkr(j+1, :);
    p4.MarkerSize = 8;
    p4.LineWidth = 1.5;
	Lgd{Lgd_cnt} = strcat( 'D gamma=', num2str(gamma_D_list(j)) );
    Lgd_cnt = Lgd_cnt + 1;
    hold on
    
    energy_max = max([ energy_line_dim(energy_line_dim~=Inf), energy_max ]);
    energy_min = min([ energy_line_dim, energy_min ]);

end

%% ------------------------------------------
% if energy_max < emergy_max_tmp
%     energy_max = emergy_max_tmp;
% end
% if energy_min > emergy_min_tmp
%     energy_min = emergy_min_tmp;
% end
hold off;
set(gca,'XTick', C_max_list);
ax = gca;
ax.YAxis.Exponent = 3;
ax.Box = 'off';
axis([min(C_max_list), max(C_max_list), 0.99*energy_min, 1.01*energy_max]);
grid on
legend(Lgd);
xlabel('Convergence Error Limit $C_{max}$','Interpreter','latex');
ylabel('Energy Cost','Interpreter','latex');


%% =========================
figure(2);
Lgd_cnt = 1;
load(strcat(pathname_data, 'Energy_opt_stepsize.mat'));
energy_min = min(energy_line_opt_stepsize);
energy_max = max(energy_line_opt_stepsize);
% ------------ Plot 1 ---------
p1 = semilogy(C_max_list, energy_line_opt_stepsize);
%         p2 = loglog(C_max_list, energy_line_opt_stepsize);
p1.Color = clr(1, :);
p1.LineStyle = str(1, :);
p1.Marker = mkr(1, :);
p1.MarkerSize = 10;
p1.LineWidth = 1.5;
Lgd{Lgd_cnt} = 'Optimal';
Lgd_cnt = Lgd_cnt + 1;
hold on;

for j = 1 : length(gamma_C_list)
    load(strcat(pathname_data, 'Energy_cons_gamma=', num2str(gamma_C_list(j)), '.mat'));
    energy_min = min([energy_min, min(energy_line_cons)]);
    energy_max = max([energy_max, max(energy_line_cons)]);
    % ------------ Plot i ---------
    p2 = semilogy(C_max_list, energy_line_cons);
    %         p1 = loglog(C_max_list, energy_line_cons);
    p2.Color = clr(2, :);
    p2.LineStyle = str(2, :);
    p2.Marker = mkr(2, :);
    p2.MarkerSize = 10;
    p2.LineWidth = 1.5;
	Lgd{Lgd_cnt} = strcat( 'Constant' );
    Lgd_cnt = Lgd_cnt + 1;
    hold on
end

for j = 1 : length(gamma_E_list)
    load(strcat(pathname_data, 'Energy_exp_gamma=', num2str(gamma_E_list(j)), '.mat'));
    energy_min = min([energy_min, min(energy_line_exp)]);
    energy_max = max([energy_max, max(energy_line_exp)]);
    p3 = semilogy(C_max_list, energy_line_exp);
    %         p1 = loglog(C_max_list, energy_line_cons);
    p3.Color = clr(3, :);
    p3.LineStyle = str(3, :);
    p3.Marker = mkr(3, :);
    p3.MarkerSize = 10;
    p3.LineWidth = 1.5;
	Lgd{Lgd_cnt} = strcat( 'Exponential' );
    Lgd_cnt = Lgd_cnt + 1;
    hold on
end

for j = 1 : length(gamma_D_list)
    load(strcat(pathname_data, 'Energy_dim_gamma=', num2str(gamma_D_list(j)), '.mat'));
    energy_min = min([energy_min, min(energy_line_dim)]);
    energy_max = max([energy_max, max(energy_line_dim)]);
    % ------------ Plot i ---------
    p4 = semilogy(C_max_list, energy_line_dim);
    %         p4 = loglog(C_max_list, energy_line_cons);
    p4.Color = clr(4, :);
    p4.LineStyle = str(4, :);
    p4.Marker = mkr(4, :);
    p4.MarkerSize = 10;
    p4.LineWidth = 1.5;
    Lgd{Lgd_cnt} = strcat( 'Diminishing' );
    Lgd_cnt = Lgd_cnt + 1;
    hold on
end

hold off;
set(gca,'XTick', C_max_list);
ax = gca;
ax.YAxis.Exponent = 3;
ax.Box = 'off';
axis([min(C_max_list), max(C_max_list), 0.99*energy_min, 1.01*energy_max]);
grid on
set(gca, 'linewidth', 1.5);
set(gca, 'FontSize', 12);
legend(Lgd);
xlabel('Convergence Error Limit $C_{max}$','Interpreter','latex');
ylabel('Energy Cost','Interpreter','latex');
