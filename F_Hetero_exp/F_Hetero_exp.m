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
pathname_tmp = 'F:\Research\MyPaper\GenQSGD_Journal\Simulation\F_Hetero_exp\Tmp\';
pathname_data = 'F:\Research\MyPaper\GenQSGD_Journal\Simulation\F_Hetero_exp\Data\';
rand_stop = 20;

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
T_max = 100000;
gamma_exp_list = [0.02];
rho_E = 0.9995;
gamma_PM_list = [0.02];
gamma_FedAvg_list = [0.02];
gamma_PR_list = [0.02];

%  ------- Setting for System Parameter ---------
F_ratio_list = [ 1:2:10 ];
F_mean = 1e9;

% ----------- Setting for Plot ---------------
figure(1);
Lgd_cnt =1;
energy_max = 0;
energy_min = 1000000;

%% ========= Simulation 1: Proposed Exponential ========
gamma_exp_list = [0.02];
rho_E = 0.9995;
for j = 1 : length(gamma_exp_list)
    gamma = gamma_exp_list(j);
    % --------- Simulation 1 ---------------
    energy_line_exp = [];
    fprintf('Simulation : Proposed Exponential\n');
    for i = 1 : length(F_ratio_list)
        F_ratio = F_ratio_list(i);
        F_min = 2*F_mean/(1+F_ratio);
        F_max = 2*F_ratio*F_mean/(1+F_ratio);
        F_n = F_max*ones(floor(N/2), 1);
        F_n = [F_n; F_min*ones(ceil(N/2), 1)];
        fprintf('F_ratio=%d start ...\n', F_ratio);
        [ energy_exp ] = E_Q_exp( N, s_n, s_0, D, alpha, ...
            C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
            Q, gamma, rho_E, G, L, sigma, C_max, T_max );
        fprintf('energy_opt=%0.1f\n',  energy_exp);
        energy_line_exp = [energy_line_exp, energy_exp];
    end
    save(strcat(pathname_data, 'Energy_exp_gamma=', num2str(gamma_exp_list(j)), '.mat'), 'energy_line_exp');
    
    % ------------ Plot 1 ---------
    p1 = plot(F_ratio_list, energy_line_exp);
    %         p2 = loglog(F_ratio_list, energy_line_opt_stepsize);
    p1.Color = clr(1, :);
    p1.LineStyle = str(1, :);
    p1.Marker = mkr(1, :);
    p1.MarkerSize = 8;
    p1.LineWidth = 1.5;
    Lgd{Lgd_cnt} = 'Proposed';
    Lgd_cnt = Lgd_cnt + 1;
    hold on
    
    energy_max = max([ energy_line_exp(energy_line_exp~=Inf) ]);
    energy_min = min([energy_line_exp]);
end

%% ========== Simulation 2: PM-SGD =============
% ---------- Setting for Simulation --------
gamma_PM_list = [0.02];

for j = 1 : length(gamma_PM_list)
    gamma_PM = gamma_PM_list(j);
    % --------- Simulation 2 ---------------
    energy_line_PM = [];
    fprintf(strcat('Simulation : PM-SGD (gamma=', num2str(gamma_PM), ')\n'));
    for i = 1 : length(F_ratio_list)
        F_ratio = F_ratio_list(i);
        F_min = 2*F_mean/(1+F_ratio);
        F_max = 2*F_ratio*F_mean/(1+F_ratio);
        F_n = F_max*ones(floor(N/2), 1);
        F_n = [F_n; F_min*ones(ceil(N/2), 1)];
        fprintf('F_ratio=%d start ...\n', F_ratio);
        [ energy_exp] = E_Q_PM( N, s_n, s_0, D, alpha, ...
            C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
            Q, gamma_PM, rho_E, G, L, sigma, C_max, T_max );
        fprintf('energyPM=%0.1f\n',  energy_exp);
        energy_line_PM = [energy_line_PM, energy_exp];
    end
    save(strcat(pathname_data, 'Energy_PM_gamma=', num2str(gamma_PM_list(j)), '.mat'), 'energy_line_PM');
    
    % ------------ Plot i ---------
    p2 = plot(F_ratio_list, energy_line_PM);
    %         p2 = loglog(F_ratio_list, energy_line_PM);
    p2.Color = clr(2, :);
    p2.LineStyle = str(2, :);
    p2.Marker = mkr(j+1, :);
    p2.MarkerSize = 8;
    p2.LineWidth = 1.5;
	Lgd{Lgd_cnt} = strcat( 'PM-SGD' );
    Lgd_cnt = Lgd_cnt + 1;
    hold on
    
    energy_max = max([ energy_line_PM(energy_line_PM~=Inf), energy_max ]);
    energy_min = min([ energy_line_PM, energy_min ]);
end

%% ========== Simulation 3: FedAvg =============
% ---------- Setting for Simulation --------
gamma_FedAvg_list = [0.02];

for j = 1 : length(gamma_FedAvg_list)
    gamma_FedAvg = gamma_FedAvg_list(j);
    % --------- Simulation 3 ---------------
    energy_line_FedAvg = [];
    fprintf(strcat('Simulation : FedAvg (gamma=', num2str(gamma_FedAvg), ')\n'));
    for i = 1 : length(F_ratio_list)
        F_ratio = F_ratio_list(i);
        F_min = 2*F_mean/(1+F_ratio);
        F_max = 2*F_ratio*F_mean/(1+F_ratio);
        F_n = F_max*ones(floor(N/2), 1);
        F_n = [F_n; F_min*ones(ceil(N/2), 1)];
        fprintf('F_ratio=%d start ...\n', F_ratio);
        [ energy_FedAvg] = E_Q_FedAvg( N, s_n, s_0, D, alpha, ...
            C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
            Q, gamma_FedAvg, rho_E, G, L, sigma, C_max, T_max );
        fprintf('energy_FedAvg=%0.1f\n',  energy_FedAvg);
        energy_line_FedAvg = [energy_line_FedAvg, energy_FedAvg];
    end
    save(strcat(pathname_data, 'Energy_FedAvg_gamma=', num2str(gamma_FedAvg_list(j)), '.mat'), 'energy_line_FedAvg');
    
    % ------------ Plot i ---------
    p3 = plot(F_ratio_list, energy_line_FedAvg);
    %         p3 = loglog(F_ratio_list, energy_line_cons);
    p3.Color = clr(3, :);
    p3.LineStyle = str(3, :);
    p3.Marker = mkr(j+1, :);
    p3.MarkerSize = 8;
    p3.LineWidth = 1.5;
	Lgd{Lgd_cnt} = strcat( 'FedAvg' );
    Lgd_cnt = Lgd_cnt + 1;
    hold on
    
    energy_max = max([ energy_line_FedAvg(energy_line_FedAvg~=Inf), energy_max ]);
    energy_min = min([ energy_line_FedAvg, energy_min ]);

end

%% ========== Simulation 4: PR-SGD =============
% ---------- Setting for Simulation --------
gamma_PR_list = [0.02];

for j = 1 : length(gamma_PR_list)
    gamma_PR = gamma_PR_list(j);
    % --------- Simulation 3 ---------------
    energy_line_PR = [];
    fprintf(strcat('Simulation : PR-SGD (gamma=', num2str(gamma_PR), ')\n'));
    for i = 1 : length(F_ratio_list)
        F_ratio = F_ratio_list(i);
        F_min = 2*F_mean/(1+F_ratio);
        F_max = 2*F_ratio*F_mean/(1+F_ratio);
        F_n = F_max*ones(floor(N/2), 1);
        F_n = [F_n; F_min*ones(ceil(N/2), 1)];
        fprintf('F_ratio=%d start ...\n', F_ratio);
        [ energy_PR ] = E_Q_PR( N, s_n, s_0, D, alpha, ...
            C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
            Q, gamma_PR, rho_E, G, L, sigma, C_max, T_max );
        fprintf('energy_PR=%0.1f\n',  energy_PR);
        energy_line_PR = [energy_line_PR, energy_PR];
    end
    save(strcat(pathname_data, 'Energy_PR_gamma=', num2str(gamma_PR_list(j)), '.mat'), 'energy_line_PR');
    
    % ------------ Plot i ---------
    p4 = plot(F_ratio_list, energy_line_PR);
    %         p4 = loglog(F_ratio_list, energy_line_cons);
    p4.Color = clr(4, :);
    p4.LineStyle = str(4, :);
    p4.Marker = mkr(j+1, :);
    p4.MarkerSize = 8;
    p4.LineWidth = 1.5;
	Lgd{Lgd_cnt} = strcat( 'PR-SGD' );
    Lgd_cnt = Lgd_cnt + 1;
    hold on
    
    energy_max = max([ energy_line_PR(energy_line_PR~=Inf), energy_max ]);
    energy_min = min([ energy_line_PR, energy_min ]);
end

%% ------------------------------------------
% if energy_max < emergy_max_tmp
%     energy_max = emergy_max_tmp;
% end
% if energy_min > emergy_min_tmp
%     energy_min = emergy_min_tmp;
% end
hold off;
set(gca,'XTick', F_ratio_list);
ax = gca;
ax.YAxis.Exponent = 3;
ax.Box = 'off';
axis([min(F_ratio_list), max(F_ratio_list), 0.99*energy_min, 1.01*energy_max]);
grid on
legend(Lgd);
xlabel('System Heterogeneity $F_{n_1}/F_{n_2}$','Interpreter','latex');
ylabel('Energy Cost','Interpreter','latex');


% %% =========================
% figure(2);
% Lgd_cnt = 1;
% for j = 1 : length(gamma_exp_list)
%     load(strcat(pathname_data, 'Energy_exp_gamma=', num2str(gamma_exp_list(j)), '.mat'));
%     energy_min = min(energy_line_exp);
%     energy_max = max(energy_line_exp);
%     % ------------ Plot 1 ---------
%     p1 = semilogy(F_ratio_list, energy_line_exp);
%     %         p2 = loglog(F_ratio_list, energy_line_cons);
%     p1.Color = clr(1, :);
%     p1.LineStyle = str(1, :);
%     p1.Marker = mkr(1, :);
%     p1.MarkerSize = 10;
%     p1.LineWidth = 1.5;
%     Lgd{Lgd_cnt} = 'GenE';
%     Lgd_cnt = Lgd_cnt + 1;
%     hold on;
% end
% 
% for j = 1 : length(gamma_PM_list)
%     load(strcat(pathname_data, 'Energy_PM_gamma=', num2str(gamma_PM_list(j)), '.mat'));
%     energy_min = min([energy_min, min(energy_line_PM)]);
%     energy_max = max([energy_max, max(energy_line_PM)]);
%     % ------------ Plot i ---------
%     p2 = semilogy(F_ratio_list, energy_line_PM);
%     %         p1 = loglog(F_ratio_list, energy_line_cons);
%     p2.Color = clr(2, :);
%     p2.LineStyle = str(2, :);
%     p2.Marker = mkr(2, :);
%     p2.MarkerSize = 10;
%     p2.LineWidth = 1.5;
% 	Lgd{Lgd_cnt} = strcat( 'PM-E' );
%     Lgd_cnt = Lgd_cnt + 1;
%     hold on
% end
% 
% for j = 1 : length(gamma_FedAvg_list)
%     load(strcat(pathname_data, 'Energy_FedAvg_gamma=', num2str(gamma_FedAvg_list(j)), '.mat'));
%     energy_min = min([energy_min, min(energy_line_FedAvg)]);
%     energy_max = max([energy_max, max(energy_line_FedAvg)]);
%     p3 = semilogy(F_ratio_list, energy_line_FedAvg);
%     %         p3 = loglog(F_ratio_list, energy_line_FedAvg);
%     p3.Color = clr(3, :);
%     p3.LineStyle = str(3, :);
%     p3.Marker = mkr(3, :);
%     p3.MarkerSize = 10;
%     p3.LineWidth = 1.5;
% 	Lgd{Lgd_cnt} = strcat( 'FA-E' );
%     Lgd_cnt = Lgd_cnt + 1;
%     hold on
% end
% 
% for j = 1 : length(gamma_PR_list)
%     load(strcat(pathname_data, 'Energy_PR_gamma=', num2str(gamma_PR_list(j)), '.mat'));
%     energy_min = min([energy_min, min(energy_line_PR)]);
%     energy_max = max([energy_max, max(energy_line_PR)]);
%     % ------------ Plot i ---------
%     p4 = semilogy(F_ratio_list, energy_line_PR);
%     %         p4 = loglog(F_ratio_list, energy_line_PR);
%     p4.Color = clr(4, :);
%     p4.LineStyle = str(4, :);
%     p4.Marker = mkr(4, :);
%     p4.MarkerSize = 10;
%     p4.LineWidth = 1.5;
%     Lgd{Lgd_cnt} = strcat( 'PR-E' );
%     Lgd_cnt = Lgd_cnt + 1;
%     hold on
% end
% 
% hold off;
% set(gca,'XTick', F_ratio_list);
% ax = gca;
% ax.YAxis.Exponent = 3;
% ax.Box = 'off';
% axis([min(F_ratio_list), max(F_ratio_list), 0.99*energy_min, 1.01*energy_max]);
% xlabel('System Heterogeneity $F_{n_1}/F_{n_2}$','Interpreter','latex');
% ylabel('Energy Cost','Interpreter','latex');
% set(gca,'linewidth',1.5);
% set(gca,'FontSize',12);
% legend(Lgd);
% grid on
%% ===============
figure(2);
load(strcat(pathname_data, 'Energy_Gen.mat'));
load(strcat(pathname_data, 'Energy_PM.mat'));
load(strcat(pathname_data, 'Energy_FedAvg.mat'));
load(strcat(pathname_data, 'Energy_PR.mat'));
line1=semilogy(F_ratio_list, energy_line_Gen,'r-o',F_ratio_list, energy_line_PM,'b:+',F_ratio_list, energy_line_FedAvg,'g:s',F_ratio_list, energy_line_PR,'m:^');
% grid on;
hold on;
set(line1,'LineWidth',1.5,'MarkerSize',8);

load(strcat(pathname_data, 'Energy_exp.mat'));
load(strcat(pathname_data, 'Energy_PM_fix.mat'));
load(strcat(pathname_data, 'Energy_FedAvg_fix.mat'));
load(strcat(pathname_data, 'Energy_PR_fix.mat'));
line2=semilogy(F_ratio_list, energy_line_exp,'r--o',F_ratio_list, energy_line_PM_fix,'b--+',F_ratio_list, energy_line_FedAvg_fix,'g--s',F_ratio_list, energy_line_PR_fix,'m--^');
set(line2,'LineWidth',1.5,'MarkerSize',8);

energy_min = min(energy_line_Gen);  energy_max = max(energy_line_Gen);
energy_min = min([energy_min, min(energy_line_exp)]);  energy_max = max([energy_max, max(energy_line_exp)]);
energy_min = min([energy_min, min(energy_line_PM)]);    energy_max = max([energy_max, max(energy_line_PM)]);
energy_min = min([energy_min, min(energy_line_PM_fix)]);    energy_max = max([energy_max, max(energy_line_PM_fix)]);
energy_min = min([energy_min, min(energy_line_FedAvg)]);    energy_max = max([energy_max, max(energy_line_FedAvg)]);
energy_min = min([energy_min, min(energy_line_FedAvg_fix)]);    energy_max = max([energy_max, max(energy_line_FedAvg_fix)]);
energy_min = min([energy_min, min(energy_line_PR)]);    energy_max = max([energy_max, max(energy_line_PR)]);
energy_min = min([energy_min, min(energy_line_PR_fix)]);    energy_max = max([energy_max, max(energy_line_PR_fix)]);

set(gca,'XTick', F_ratio_list);
ax = gca;
ax.YAxis.Exponent = 5;
ax.Box = 'off';
axis([min(F_ratio_list), max(F_ratio_list), 0.99*energy_min, 1.5*energy_max]);
xlabel('Difference of Computation Capabilities $F^{(1)}/F^{(2)}$','Interpreter','latex');
ylabel('Energy Cost','Interpreter','latex');
set(gca,'linewidth',1.5);
set(gca,'FontSize',12);
grid on;

lgd1 = legend(line1,'Gen-O' , 'PM-E-opt' , 'FA-E-opt' , 'PR-E-opt', 'Without Reflector','Location','northwest');
ah=axes('position',get(gca,'position'),'visible','off');
lgd2 = legend(ah,line2,'Gen-E' ,  'PM-E-fix' , 'FA-E-fix' , 'PR-E-fix', 'Location','northeast');
set(lgd1,'FontSize',10, 'LineWidth',1);
set(lgd2,'FontSize',10, 'LineWidth',1);