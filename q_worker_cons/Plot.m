pathname_tmp = '.\Tmp\';
pathname_data = '.\Data\';

bits_list = [10:16, 20, 24, 28, 32];

figure(1);
load(strcat(pathname_data, 'Energy_Gen.mat'));
load(strcat(pathname_data, 'Energy_PM.mat'));
load(strcat(pathname_data, 'Energy_FedAvg.mat'));
load(strcat(pathname_data, 'Energy_PR.mat'));
% line1=semilogy(bits_list, energy_line_Gen,'r-o',bits_list, energy_line_PM,'b:+',bits_list, energy_line_FedAvg,'g:s',bits_list, energy_line_PR,'m:^');
line1=loglog(bits_list(2:end), energy_line_Gen(2:end),'r-o',bits_list(2:end), energy_line_PM(2:end),'b:+',bits_list(2:end), energy_line_FedAvg(2:end),'g:s',bits_list(2:end), energy_line_PR(2:end),'m:^');

% grid on;
hold on;
set(line1,'LineWidth',1.5,'MarkerSize',8);

load(strcat(pathname_data, 'Energy_cons.mat'));
load(strcat(pathname_data, 'Energy_PM_fix.mat'));
load(strcat(pathname_data, 'Energy_FedAvg_fix.mat'));
load(strcat(pathname_data, 'Energy_PR_fix.mat'));
line2=semilogy(bits_list(2:end), energy_line_cons(2:end),'r--o',bits_list(2:end), energy_line_PM_fix(2:end),'b--+',bits_list(2:end), energy_line_FedAvg_fix(2:end),'g--s',bits_list(2:end), energy_line_PR_fix(2:end),'m--^');
set(line2,'LineWidth',1.5,'MarkerSize',8);

energy_min = min(energy_line_Gen);  energy_max = max(energy_line_Gen);
energy_min = min([energy_min, min(energy_line_cons)]);  energy_max = max([energy_max, max(energy_line_cons)]);
energy_min = min([energy_min, min(energy_line_PM)]);    energy_max = max([energy_max, max(energy_line_PM)]);
energy_min = min([energy_min, min(energy_line_PM_fix)]);    energy_max = max([energy_max, max(energy_line_PM_fix)]);
energy_min = min([energy_min, min(energy_line_FedAvg)]);    energy_max = max([energy_max, max(energy_line_FedAvg)]);
energy_min = min([energy_min, min(energy_line_FedAvg_fix)]);    energy_max = max([energy_max, max(energy_line_FedAvg_fix)]);
energy_min = min([energy_min, min(energy_line_PR)]);    energy_max = max([energy_max, max(energy_line_PR)]);
energy_min = min([energy_min, min(energy_line_PR_fix)]);    energy_max = max([energy_max, max(energy_line_PR_fix)]);

set(gca,'XTick', bits_list);
ax = gca;
ax.YAxis.Exponent = 4;
ax.Box = 'off';
axis([min(bits_list(2:end)), max(bits_list), 0.99*energy_min, 1.5*energy_max]);
xlabel('Quantization Parameters $\log_2s_n,n\in\mathcal N$','Interpreter','latex');
ylabel('Energy Cost','Interpreter','latex');
set(gca,'linewidth',1.5);
set(gca,'FontSize',12);
grid on;

lgd1 = legend(line1,'Gen-O' , 'PM-C-opt' , 'FA-C-opt' , 'PR-C-opt', 'Without Reflector','Location','northwest');
ah=axes('position',get(gca,'position'),'visible','off');
lgd2 = legend(ah,line2,'Gen-C' ,  'PM-C-fix' , 'FA-C-fix' , 'PR-C-fix', 'Location','northeast');
set(lgd1,'FontSize',10, 'LineWidth',1);
set(lgd2,'FontSize',10, 'LineWidth',1);