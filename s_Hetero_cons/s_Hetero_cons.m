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
bits_0 = [32];
s_0 = 2^bits_0;

D = 128*784;
Q = 2.5;          % 3.2-0.7=2.5
G = 33.63;      % 0.6
L = 0.084;     %1.4e3~1.5e3        %2.4163     22   L=1
sigma = 33.18;     % 18

clr=[[1 0 0];[1 0 1];[0 0 1];[0 0.5 0.5];[0.2 0.5 0.3];[0 0 0]];
str=[' -';' :';'--'; '-.'];
mkr = ['o';'+';'s';'^';'*';'d'];
pathname_tmp = 'F:\Research\MyPaper\GenQSGD_Journal\Simulation\s_Hetero_cons\Tmp\';
pathname_data = 'F:\Research\MyPaper\GenQSGD_Journal\Simulation\s_Hetero_cons\Data\';
rand_stop = 20;

%% ============ Setting for Step Size Rules: Varying Computation Heterogeneous =============
%  ------- Setting for System Parameter ---------
F_max = 1e9;
F_n = F_max*ones(N, 1);
C_n_max = 100e6;    % 100e6
C_n = C_n_max*ones(N, 1);

bits_n = [14];
s_mean = 2^bits_n;
s_ratio_list = 1:2:10;

p_mean = 1.5;
r_mean = 50e6 / N;
p_n = p_mean*ones(N,1);
r_n = r_mean*ones(N,1);

% ---------- Setting for Simulation --------
C_max = 0.25;
T_max = 100000;
gamma_cons_list = [0.01];
gamma_PM_list = [0.01];
gamma_FedAvg_list = [0.01];
gamma_PR_list = [0.01];
% ----------- Setting for Plot ---------------
figure(1);
Lgd_cnt =1;
energy_max = 0;
energy_min = 100000;

%% ========= Simulation 1: Proposed Constant ========
gamma_cons_list = [0.01];
for j = 1 : length(gamma_cons_list)
    gamma = gamma_cons_list(j);
    % --------- Simulation 1 ---------------
    energy_line_cons = [];
    fprintf('Simulation : Proposed Constant\n');
    for i = 1 : length(s_ratio_list)
        s_ratio = s_ratio_list(i);
        s_min = 2*s_mean/(1+s_ratio);
        s_max = 2*s_ratio*s_mean/(1+s_ratio);
        s_n = s_max*ones(floor(N/2), 1);
        s_n = [s_n; s_min*ones(ceil(N/2), 1)];
        fprintf('s_ratio=%d start ...\n', s_ratio);
        [ energy_cons ] = E_Q_cons( N, s_n, s_0, D, alpha, ...
            C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
            Q, gamma, G, L, sigma, C_max, T_max );
        fprintf('energy_opt=%0.1f\n',  energy_cons);
        energy_line_cons = [energy_line_cons, energy_cons];
    end
    save(strcat(pathname_data, 'Energy_cons_gamma=', num2str(gamma_cons_list(j)), '.mat'), 'energy_line_cons');
    
    % ------------ Plot 1 ---------
    p1 = plot(s_ratio_list, energy_line_cons);
    %         p2 = loglog(s_ratio_list, energy_line_opt_stepsize);
    p1.Color = clr(1, :);
    p1.LineStyle = str(1, :);
    p1.Marker = mkr(1, :);
    p1.MarkerSize = 8;
    p1.LineWidth = 1.5;
    Lgd{Lgd_cnt} = 'Proposed';
    Lgd_cnt = Lgd_cnt + 1;
    hold on
    
    energy_max = max([ energy_line_cons(energy_line_cons~=Inf) ]);
    energy_min = min([energy_line_cons]);
end

%% ========== Simulation 2: PM-SGD =============
% ---------- Setting for Simulation --------
gamma_PM_list = [0.01];
for j = 1 : length(gamma_PM_list)
    gamma_PM = gamma_PM_list(j);
    % --------- Simulation 2 ---------------
    energy_line_PM = [];
    fprintf(strcat('Simulation : PM-SGD (gamma=', num2str(gamma_PM), ')\n'));
    for i = 1 : length(s_ratio_list)
        s_ratio = s_ratio_list(i);
        s_min = 2*s_mean/(1+s_ratio);
        s_max = 2*s_ratio*s_mean/(1+s_ratio);
        s_n = s_max*ones(floor(N/2), 1);
        s_n = [s_n; s_min*ones(ceil(N/2), 1)];
        fprintf('s_ratio=%d start ...\n', s_ratio);
        [ energy_cons] = E_Q_PM( N, s_n, s_0, D, alpha, ...
            C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
            Q, gamma_PM, G, L, sigma, C_max, T_max );
        fprintf('energy_cons=%0.1f\n',  energy_cons);
        energy_line_PM = [energy_line_PM, energy_cons];
    end
    save(strcat(pathname_data, 'Energy_PM_gamma=', num2str(gamma_PM_list(j)), '.mat'), 'energy_line_PM');
    
    % ------------ Plot i ---------
    p2 = plot(s_ratio_list, energy_line_PM);
    %         p2 = loglog(s_ratio_list, energy_line_PM);
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
gamma_FedAvg_list = [0.01];

for j = 1 : length(gamma_FedAvg_list)
    gamma_FedAvg = gamma_FedAvg_list(j);
    % --------- Simulation 3 ---------------
    energy_line_FedAvg = [];
    fprintf(strcat('Simulation : FedAvg (gamma=', num2str(gamma_FedAvg), ')\n'));
    for i = 1 : length(s_ratio_list)
        s_ratio = s_ratio_list(i);
        s_min = 2*s_mean/(1+s_ratio);
        s_max = 2*s_ratio*s_mean/(1+s_ratio);
        s_n = s_max*ones(floor(N/2), 1);
        s_n = [s_n; s_min*ones(ceil(N/2), 1)];
        fprintf('s_ratio=%d start ...\n', s_ratio);
        [ energy_FedAvg] = E_Q_FedAvg( N, s_n, s_0, D, alpha, ...
            C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
            Q, gamma_FedAvg, G, L, sigma, C_max, T_max );
        fprintf('energy_FedAvg=%0.1f\n',  energy_FedAvg);
        energy_line_FedAvg = [energy_line_FedAvg, energy_FedAvg];
    end
    save(strcat(pathname_data, 'Energy_FedAvg_gamma=', num2str(gamma_FedAvg_list(j)), '.mat'), 'energy_line_FedAvg');
    
    % ------------ Plot i ---------
    p3 = plot(s_ratio_list, energy_line_FedAvg);
    %         p3 = loglog(s_ratio_list, energy_line_cons);
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
gamma_PR_list = [0.01];

for j = 1 : length(gamma_PR_list)
    gamma_PR = gamma_PR_list(j);
    % --------- Simulation 3 ---------------
    energy_line_PR = [];
    fprintf(strcat('Simulation : PR-SGD (gamma=', num2str(gamma_PR), ')\n'));
    for i = 1 : length(s_ratio_list)
        s_ratio = s_ratio_list(i);
        s_min = 2*s_mean/(1+s_ratio);
        s_max = 2*s_ratio*s_mean/(1+s_ratio);
        s_n = s_max*ones(floor(N/2), 1);
        s_n = [s_n; s_min*ones(ceil(N/2), 1)];
        fprintf('s_ratio=%d start ...\n', s_ratio);
        [ energy_PR ] = E_Q_PR( N, s_n, s_0, D, alpha, ...
            C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
            Q, gamma_PR, G, L, sigma, C_max, T_max );
        fprintf('energy_PR=%0.1f\n',  energy_PR);
        energy_line_PR = [energy_line_PR, energy_PR];
    end
    save(strcat(pathname_data, 'Energy_PR_gamma=', num2str(gamma_PR_list(j)), '.mat'), 'energy_line_PR');
    
    % ------------ Plot i ---------
    p4 = plot(s_ratio_list, energy_line_PR);
    %         p4 = loglog(s_ratio_list, energy_line_cons);
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
set(gca,'XTick', s_ratio_list);
ax = gca;
ax.YAxis.Exponent = 3;
ax.Box = 'off';
axis([min(s_ratio_list), max(s_ratio_list), 0.99*energy_min, 1.01*energy_max]);
grid on
legend(Lgd);
xlabel('System Heterogeneity $s_{n_1}/s_{n_2}$','Interpreter','latex');
ylabel('Energy Cost','Interpreter','latex');


% %% =========================
% figure(2);
% Lgd_cnt = 1;
% for j = 1 : length(gamma_cons_list)
%     load(strcat(pathname_data, 'Energy_cons_gamma=', num2str(gamma_cons_list(j)), '.mat'));
%     energy_min = min(energy_line_cons);
%     energy_max = max(energy_line_cons);
%     % ------------ Plot 1 ---------
% %     p1 = semilogy(s_ratio_list, energy_line_cons);
%     p1 = plot(s_ratio_list, energy_line_cons);
% %             p1 = loglog(s_ratio_list, energy_line_cons);
%     p1.Color = clr(1, :);
%     p1.LineStyle = str(1, :);
%     p1.Marker = mkr(1, :);
%     p1.MarkerSize = 10;
%     p1.LineWidth = 1.5;
%     Lgd{Lgd_cnt} = 'Gen-C';
%     Lgd_cnt = Lgd_cnt + 1;
%     hold on;
% end
% 
% for j = 1 : length(gamma_PM_list)
%     load(strcat(pathname_data, 'Energy_PM_gamma=', num2str(gamma_PM_list(j)), '.mat'));
%     energy_min = min([energy_min, min(energy_line_PM)]);
%     energy_max = max([energy_max, max(energy_line_PM)]);
%     % ------------ Plot i ---------
% %     p2 = semilogy(s_ratio_list, energy_line_PM);
%     p2 = plot(s_ratio_list, energy_line_PM);
% %             p2 = loglog(s_ratio_list, energy_line_cons);
%     p2.Color = clr(2, :);
%     p2.LineStyle = str(2, :);
%     p2.Marker = mkr(2, :);
%     p2.MarkerSize = 10;
%     p2.LineWidth = 1.5;
% 	Lgd{Lgd_cnt} = strcat( 'PM-C' );
%     Lgd_cnt = Lgd_cnt + 1;
%     hold on
% end
% 
% for j = 1 : length(gamma_FedAvg_list)
%     load(strcat(pathname_data, 'Energy_FedAvg_gamma=', num2str(gamma_FedAvg_list(j)), '.mat'));
%     energy_min = min([energy_min, min(energy_line_FedAvg)]);
%     energy_max = max([energy_max, max(energy_line_FedAvg)]);
% %     p3 = semilogy(s_ratio_list, energy_line_FedAvg);
%     p3 = plot(s_ratio_list, energy_line_FedAvg);
% %             p3 = loglog(s_ratio_list, energy_line_FedAvg);
%     p3.Color = clr(3, :);
%     p3.LineStyle = str(3, :);
%     p3.Marker = mkr(3, :);
%     p3.MarkerSize = 10;
%     p3.LineWidth = 1.5;
% 	Lgd{Lgd_cnt} = strcat( 'FA-C' );
%     Lgd_cnt = Lgd_cnt + 1;
%     hold on
% end
% 
% for j = 1 : length(gamma_PR_list)
%     load(strcat(pathname_data, 'Energy_PR_gamma=', num2str(gamma_PR_list(j)), '.mat'));
%     energy_min = min([energy_min, min(energy_line_PR)]);
%     energy_max = max([energy_max, max(energy_line_PR)]);
%     % ------------ Plot i ---------
% %     p4 = semilogy(s_ratio_list, energy_line_PR);
%     p4 = plot(s_ratio_list, energy_line_PR);
% %             p4 = loglog(s_ratio_list, energy_line_PR);
%     p4.Color = clr(4, :);
%     p4.LineStyle = str(4, :);
%     p4.Marker = mkr(4, :);
%     p4.MarkerSize = 10;
%     p4.LineWidth = 1.5;
%     Lgd{Lgd_cnt} = strcat( 'PR-C' );
%     Lgd_cnt = Lgd_cnt + 1;
%     hold on
% end
% 
% hold off;
% set(gca,'XTick', s_ratio_list);
% ax = gca;
% ax.YAxis.Exponent = 3;
% ax.Box = 'off';
% axis([min(s_ratio_list), max(s_ratio_list), 0.99*energy_min, 1.01*energy_max]);
% xlabel('System Heterogeneity $s_{n_1}/s_{n_2}$','Interpreter','latex');
% ylabel('Energy Cost','Interpreter','latex');
% set(gca,'linewidth',1.5);
% set(gca,'FontSize',12);
% legend(Lgd);
% grid on
%% =========================
figure(2);
load(strcat(pathname_data, 'Energy_Gen.mat'));
load(strcat(pathname_data, 'Energy_PM.mat'));
load(strcat(pathname_data, 'Energy_FedAvg.mat'));
load(strcat(pathname_data, 'Energy_PR.mat'));
line1=semilogy(s_ratio_list, energy_line_Gen,'r-o',s_ratio_list, energy_line_PM,'b:+',s_ratio_list, energy_line_FedAvg,'g:s',s_ratio_list, energy_line_PR,'m:^');
% grid on;
hold on;
set(line1,'LineWidth',1.5,'MarkerSize',8);

load(strcat(pathname_data, 'Energy_cons.mat'));
load(strcat(pathname_data, 'Energy_PM_fix.mat'));
load(strcat(pathname_data, 'Energy_FedAvg_fix.mat'));
load(strcat(pathname_data, 'Energy_PR_fix.mat'));
line2=semilogy(s_ratio_list, energy_line_cons,'r--o',s_ratio_list, energy_line_PM_fix,'b--+',s_ratio_list, energy_line_FedAvg_fix,'g--s',s_ratio_list, energy_line_PR_fix,'m--^');
set(line2,'LineWidth',1.5,'MarkerSize',8);

energy_min = min(energy_line_Gen);  energy_max = max(energy_line_Gen);
energy_min = min([energy_min, min(energy_line_cons)]);  energy_max = max([energy_max, max(energy_line_cons)]);
energy_min = min([energy_min, min(energy_line_PM)]);    energy_max = max([energy_max, max(energy_line_PM)]);
energy_min = min([energy_min, min(energy_line_PM_fix)]);    energy_max = max([energy_max, max(energy_line_PM_fix)]);
energy_min = min([energy_min, min(energy_line_FedAvg)]);    energy_max = max([energy_max, max(energy_line_FedAvg)]);
energy_min = min([energy_min, min(energy_line_FedAvg_fix)]);    energy_max = max([energy_max, max(energy_line_FedAvg_fix)]);
energy_min = min([energy_min, min(energy_line_PR)]);    energy_max = max([energy_max, max(energy_line_PR)]);
energy_min = min([energy_min, min(energy_line_PR_fix)]);    energy_max = max([energy_max, max(energy_line_PR_fix)]);

set(gca,'XTick', s_ratio_list);
ax = gca;
ax.YAxis.Exponent = 5;
ax.Box = 'off';
axis([min(s_ratio_list), max(s_ratio_list), 0.99*energy_min, 1.8*energy_max]);
xlabel('Difference of Quantization Parameters $s^{(1)}/s^{(2)}$','Interpreter','latex');
ylabel('Energy Cost','Interpreter','latex');
set(gca,'linewidth',1.5);
set(gca,'FontSize',12);
grid on;

lgd1 = legend(line1,'Gen-O' , 'PM-C-opt' , 'FA-C-opt' , 'PR-C-opt', 'Without Reflector','Location','northwest');
ah=axes('position',get(gca,'position'),'visible','off');
lgd2 = legend(ah,line2,'Gen-C' ,  'PM-C-fix' , 'FA-C-fix' , 'PR-C-fix', 'Location','northeast');
set(lgd1,'FontSize',10, 'LineWidth',1);
set(lgd2,'FontSize',10, 'LineWidth',1);