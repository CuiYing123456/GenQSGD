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
C_max = [0.25];
bits_list = [10:16, 20, 24, 28, 32];

% ----------- Setting for Plot ---------------
figure(1);
Lgd_cnt =1;
energy_max = 0;
energy_min = 1000000;

%% ========= Simulation 1: Proposed Exponential ========
% --------- Simulation 1 ---------------
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
save(strcat(pathname_data, 'Energy_Gen.mat'), 'energy_line_Gen');
%% =========================
figure(2);
Lgd_cnt = 1;
load(strcat(pathname_data, 'Energy_Gen.mat'));
energy_min = min(energy_line_Gen);
energy_max = max(energy_line_Gen);
% ------------ Plot 1 ---------
p1 = semilogy(bits_list, energy_line_Gen);
%         p2 = loglog(F_ratio_list, energy_line_Gen);
p1.Color = clr(1, :);
p1.LineStyle = str(1, :);
p1.Marker = mkr(1, :);
p1.MarkerSize = 10;
p1.LineWidth = 1.5;
Lgd{Lgd_cnt} = 'GenQSGD';
Lgd_cnt = Lgd_cnt + 1;
hold on;

hold off;
set(gca,'XTick', F_ratio_list);
ax = gca;
ax.YAxis.Exponent = 3;
ax.Box = 'off';
axis([min(F_ratio_list), max(F_ratio_list), 0.99*energy_min, 1.01*energy_max]);
xlabel('System Heterogeneity $\max_{n\in\mathcal N}F_n/\min_{n\in\mathcal N}F_n$','Interpreter','latex');
ylabel('Energy Cost','Interpreter','latex');
set(gca,'linewidth',1.5);
set(gca,'FontSize',12);
legend(Lgd);
grid on
