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
% F_n = 2e9+2e8*(rand(N,1)-0.5);      % 2e9
% C_n = 1000e5+100e5*(rand(N,1)-0.5);     %256*784
% p_n = 2+0.2*(rand(N,1)-0.5);
% r_n = 5e5+5e4*(rand(N,1)-0.5);

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

%% ============ Setting for Step Size Rules: Varying Computation Heterogeneous =============
F_mean = 1e9;
F_ratio = 10;
F_min = 2*F_mean/(1+F_ratio);
F_max = 2*F_ratio*F_mean/(1+F_ratio);
F_n = F_max*ones(floor(N/2), 1);
F_n = [F_n; F_min*ones(ceil(N/2), 1)];
C_n_max = 100e6;    % 100e6
C_n = C_n_max*ones(N, 1);
s_0 = 2^32;

p_mean = 1.5;
r_mean = 50e6 / N;
p_n = p_mean*ones(N,1);
r_n = r_mean*ones(N,1);

% ---------- Setting for Simulation --------
T_max = 100000;
C_max = [0.25];
bits_list = [10:16, 20, 24, 28, 32];
gamma_cons = [0.01];   %  0.02 OK!
gamma_PM = [0.01]; 
IsOpt_PM = 0;
gamma_FedAvg = [0.01];
IsOpt_FedAvg = 0;
gamma_PR = [0.01];
IsOpt_PR = 0;

energy_max = 0;
energy_min = 1000000000000;

%% ========== Simulation 2: PM-SGD =============
% --------- Simulation 2 ---------------
energy_line_PM = [];
fprintf(strcat('Simulation : PM-SGD-fix\n'));
for i = 1 : length(bits_list)
    s_n = 2^bits_list(i) * ones(N,1);
    fprintf('s_n=%d start ...\n', s_n(1));
    [ energy_PM] = E_Q_PM_fix( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_PM, G, L, sigma, C_max, T_max );
    fprintf('energy_cons=%0.1f\n',  energy_PM);
    energy_line_PM = [energy_line_PM, energy_PM];
end
energy_line_PM_fix = energy_line_PM;
save(strcat(pathname_data, 'Energy_PM_fix.mat'), 'energy_line_PM_fix');

energy_max = max([ energy_line_PM(energy_line_PM~=Inf), energy_max ]);
energy_min = min([ energy_line_PM ], energy_min);


%% ========== Simulation 3: FedAvg =============
% --------- Simulation 3 ---------------
energy_line_FedAvg = [];
fprintf(strcat('Simulation : FedAvg-fix\n'));
for i = 1 : length(bits_list)
    s_n = 2^bits_list(i) * ones(N,1);
    fprintf('s_n=%d start ...\n', s_n(1));
    [ energy_FedAvg] = E_Q_FedAvg_fix( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_FedAvg, G, L, sigma, C_max, T_max );
    fprintf('energy_FedAvg=%0.1f\n',  energy_FedAvg);
    energy_line_FedAvg = [energy_line_FedAvg, energy_FedAvg];
end
energy_line_FedAvg_fix = energy_line_FedAvg;
save(strcat(pathname_data, 'Energy_FedAvg_fix.mat'), 'energy_line_FedAvg_fix');

energy_max = max([ energy_line_FedAvg(energy_line_FedAvg~=Inf), energy_max ]);
energy_min = min([ energy_line_FedAvg, energy_min ]);


%% ========== Simulation 4: PR-SGD =============
% --------- Simulation 3 ---------------
energy_line_PR = [];
fprintf(strcat('Simulation : PR-SGD-fix\n'));
for i = 1 : length(bits_list)
    s_n = 2^bits_list(i) * ones(N,1);
    fprintf('s_n=%d start ...\n', s_n(1));
    [ energy_PR ] = E_Q_PR_fix( N, s_n, s_0, D, alpha, ...
        C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
        Q, gamma_PR, G, L, sigma, C_max, T_max );
    fprintf('energy_PR=%0.1f\n',  energy_PR);
    energy_line_PR = [energy_line_PR, energy_PR];
end
energy_line_PR_fix = energy_line_PR;
save(strcat(pathname_data, 'Energy_PR_fix.mat'), 'energy_line_PR_fix');

energy_max = max([ energy_line_PR(energy_line_PR~=Inf), energy_max ]);
energy_min = min([ energy_line_PR, energy_min ]);


%% =========================
figure(2);
Lgd_cnt = 1;
load(strcat(pathname_data, 'Energy_cons.mat'));
energy_min = min(energy_line_cons);
energy_max = max(energy_line_cons);
% ------------ Plot 1 ---------
p1 = semilogy(bits_list, energy_line_cons);
%         p2 = loglog(F_ratio_list, energy_line_cons);
p1.Color = clr(1, :);
p1.LineStyle = str(1, :);
p1.Marker = mkr(1, :);
p1.MarkerSize = 10;
p1.LineWidth = 1.5;
Lgd{Lgd_cnt} = 'Proposed';
Lgd_cnt = Lgd_cnt + 1;
hold on;


load(strcat(pathname_data, 'Energy_PM.mat'));
energy_min = min([energy_min, min(energy_line_PM)]);
energy_max = max([energy_max, max(energy_line_PM)]);
% ------------ Plot i ---------
p2 = semilogy(bits_list, energy_line_PM);
%         p1 = loglog(F_ratio_list, energy_line_cons);
p2.Color = clr(2, :);
p2.LineStyle = str(1, :);
p2.Marker = mkr(2, :);
p2.MarkerSize = 10;
p2.LineWidth = 1.5;
Lgd{Lgd_cnt} = strcat( 'PM-SGD-opt' );
Lgd_cnt = Lgd_cnt + 1;
hold on
load(strcat(pathname_data, 'Energy_PM_fix.mat'));
energy_min = min([energy_min, min(energy_line_PM_fix)]);
energy_max = max([energy_max, max(energy_line_PM_fix)]);

% ------------ Plot i ---------
p2 = semilogy(bits_list, energy_line_PM_fix);
%         p1 = loglog(F_ratio_list, energy_line_cons);
p2.Color = clr(2, :);
p2.LineStyle = str(2, :);
p2.Marker = mkr(2, :);
p2.MarkerSize = 10;
p2.LineWidth = 1.5;
Lgd{Lgd_cnt} = strcat( 'PM-SGD-fix' );
Lgd_cnt = Lgd_cnt + 1;
hold on

load(strcat(pathname_data, 'Energy_FedAvg.mat'));
energy_min = min([energy_min, min(energy_line_FedAvg)]);
energy_max = max([energy_max, max(energy_line_FedAvg)]);
p3 = semilogy(bits_list, energy_line_FedAvg);
%         p3 = loglog(F_ratio_list, energy_line_FedAvg);
p3.Color = clr(3, :);
p3.LineStyle = str(1, :);
p3.Marker = mkr(3, :);
p3.MarkerSize = 10;
p3.LineWidth = 1.5;
Lgd{Lgd_cnt} = strcat( 'FedAvg-opt' );
Lgd_cnt = Lgd_cnt + 1;
hold on
load(strcat(pathname_data, 'Energy_FedAvg_fix.mat'));
energy_min = min([energy_min, min(energy_line_FedAvg_fix)]);
energy_max = max([energy_max, max(energy_line_FedAvg_fix)]);
p3 = semilogy(bits_list, energy_line_FedAvg_fix);
%         p3 = loglog(F_ratio_list, energy_line_FedAvg);
p3.Color = clr(3, :);
p3.LineStyle = str(2, :);
p3.Marker = mkr(3, :);
p3.MarkerSize = 10;
p3.LineWidth = 1.5;
Lgd{Lgd_cnt} = strcat( 'FedAvg-fix' );
Lgd_cnt = Lgd_cnt + 1;
hold on


load(strcat(pathname_data, 'Energy_PR.mat'));
energy_min = min([energy_min, min(energy_line_PR)]);
energy_max = max([energy_max, max(energy_line_PR)]);
% ------------ Plot i ---------
p4 = semilogy(bits_list, energy_line_PR);
%         p4 = loglog(F_ratio_list, energy_line_PR);
p4.Color = clr(4, :);
p4.LineStyle = str(1, :);
p4.Marker = mkr(4, :);
p4.MarkerSize = 10;
p4.LineWidth = 1.5;
Lgd{Lgd_cnt} = strcat( 'PR-SGD-opt' );
Lgd_cnt = Lgd_cnt + 1;
load(strcat(pathname_data, 'Energy_PR_fix.mat'));
energy_min = min([energy_min, min(energy_line_PR_fix)]);
energy_max = max([energy_max, max(energy_line_PR_fix)]);
% ------------ Plot i ---------
p4 = semilogy(bits_list, energy_line_PR_fix);
%         p4 = loglog(F_ratio_list, energy_line_PR);
p4.Color = clr(4, :);
p4.LineStyle = str(2, :);
p4.Marker = mkr(4, :);
p4.MarkerSize = 10;
p4.LineWidth = 1.5;
Lgd{Lgd_cnt} = strcat( 'PR-SGD-fix' );
Lgd_cnt = Lgd_cnt + 1;


hold off;
set(gca,'XTick', bits_list);
ax = gca;
ax.YAxis.Exponent = 3;
ax.Box = 'off';
axis([min(bits_list), max(bits_list), 0.99*energy_min, 1.01*energy_max]);
xlabel('Quantization parameters $\log_2s_n,n\in\mathcal N$','Interpreter','latex');
ylabel('Energy Cost','Interpreter','latex');
set(gca,'linewidth',1.5);
set(gca,'FontSize',12);
legend(Lgd);
grid on