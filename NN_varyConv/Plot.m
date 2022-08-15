
clr=[[1 0 0];[1 0 1];[0 0 1];[0 0.5 0.5];[0.2 0.5 0.3];[0 0 0]];
str=[' -';' :';'--'; '-.'];
mkr = ['o';'+';'s';'^';'*';'d'];
pathname_data = './Data/';

G_th_list = 0.2 : 0.01 : 0.25;

figure(1);
yyaxis left;        % N-left-time
line1=plot(G_th_list, Acc_rcd_opt,'r-o',G_th_list, Acc_rcd_cons,'r-+',G_th_list, Acc_rcd_exp,'r-s',G_th_list, Acc_rcd_dim,'r-^');
grid on;
hold on;
axis([min(G_th_list), max(G_th_list), 0.95*Acc_min 1.04*Acc_max]);
ylabel('Test Accuracy ($\%$)','interpreter','latex','fontsize',15);
set(gca,'ycolor','r');
set(gca,'linewidth',1.2);
set(gca,'FontSize',15);

yyaxis right;        % N-left-time
line2=plot(G_th_list, Loss_rcd_opt,'b--o',G_th_list, Loss_rcd_cons,'b--+',G_th_list, Loss_rcd_exp,'b--s',G_th_list, Loss_rcd_dim,'b--^');
grid on;
axis([ min(G_th_list), max(G_th_list), 0.9*Loss_min 1.8*Loss_max]);
ax = gca;
ax.Box = 'off';
ylabel('Training Loss','interpreter','latex','fontsize',15);
set(gca,'ycolor','b');
set(gca,'linewidth',1.2);
set(gca,'FontSize',15);
xlabel('Convergence Error Limit $C_{\max}$','interpreter','latex','fontsize',15);

set(line1,'LineWidth',1.5,'MarkerSize',10);
set(line2,'LineWidth',1.5,'MarkerSize',10);
lgd1 = legend(line1,'accuracy (Gen-O)' , 'accuracy (Gen-C)' , 'accuracy (Gen-E)' , 'accuracy (Gen-D)', 'Without Reflector','Location','northwest');
ah=axes('position',get(gca,'position'),'visible','off');
lgd2 = legend(ah,line2,'loss (Gen-O)' , 'loss (Gen-C)' , 'loss (Gen-E)' , 'loss (Gen-D)' , 'Location','northeast');
set(lgd1,'FontSize',10, 'LineWidth',1.2);
set(lgd2,'FontSize',10, 'LineWidth',1.2);

hold off;