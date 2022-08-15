% clr = [[1, 0, 0]; [1, 0.5, 0]; [1, 1, 0]; 
%     [0, 0.5, 0]; [0.5, 1, 0]; [0.5, 1, 0.5];
%     [0, 0.2, 1]; [0, 0.8, 1]; [0.5, 1, 1]];
% clr = [[1, 0, 0]; [1, 0.5, 0]; [1, 1, 0]; 
%     [0, 0.8, 0]; [0, 1, 0.5]; [0, 1, 1];
%     [0.3, 0.2, 1]; [0.8, 0.3, 1]; [1, 0.5, 1]];
clr = [[1, 0, 0]; [1, 0.5, 0]; [1, 0.8, 0.2]; 
    [0, 0.5, 0]; [0.2, 0.8, 0]; [0.6, 1, 0];
    [0, 0.2, 1]; [0, 0.6, 1]; [0, 1, 1]];
str = {'-'; '--'; '-.'; ':'};

%% =========== Plot Test Accuracy ===========
Lgd1 = {};
figure(1);

load(sprintf('./Data/acc_opt.mat'));
p = plot(a_mean);  p.LineWidth = 1.5;  p.Color = clr(1,:);  p.LineStyle = str{1};  p.MarkerSize = 8;  hold on;
Lgd1{1} = sprintf('Gen-O');
load(sprintf('./Data/acc_cons.mat'));
p = plot(a_mean);  p.LineWidth = 1.5;  p.Color = clr(2,:);  p.LineStyle = str{1};  p.MarkerSize = 8;  hold on;
Lgd1{2} = sprintf('Gen-C');
load(sprintf('./Data/acc_exp.mat'));
p = plot(a_mean);  p.LineWidth = 1.5;  p.Color = clr(5,:);  p.LineStyle = str{1};  p.MarkerSize = 8;  hold on;
Lgd1{3} = sprintf('Gen-E');
load(sprintf('./Data/acc_dim.mat'));
p = plot(a_mean);  p.LineWidth = 1.5;  p.Color = clr(8,:);  p.LineStyle = str{1};  p.MarkerSize = 8;  hold on;
Lgd1{4} = sprintf('Gen-D');

hold off;
legend(Lgd1);
grid on;
axis([-Inf, Inf, 70, 95]);
xlabel('Global Iteration','Interpreter','latex');
ylabel('Test Accuracy (\%)','Interpreter','latex');
set(gca,'linewidth',1.5);
set(gca,'fontsize',15);
set(gca,'Box','off');

%% =========== Plot Cost Function ===========
Lgd2 = {};
figure(2);

load(sprintf('./Data/cost_opt.mat'));
p = plot(c_mean);  p.LineWidth = 1.5;  p.Color = clr(1,:);  p.LineStyle = str{1};  p.MarkerSize = 8;  hold on;
Lgd2{1} = sprintf('Gen-O');
load(sprintf('./Data/cost_cons.mat'));
p = plot(c_mean);  p.LineWidth = 1.5;  p.Color = clr(2,:);  p.LineStyle = str{1};  p.MarkerSize = 8;  hold on;
Lgd2{2} = sprintf('Gen-C');
load(sprintf('./Data/cost_exp.mat'));
p = plot(c_mean);  p.LineWidth = 1.5;  p.Color = clr(5,:);  p.LineStyle = str{1};  p.MarkerSize = 8;  hold on;
Lgd2{3} = sprintf('Gen-E');
load(sprintf('./Data/cost_dim.mat'));
p = plot(c_mean);  p.LineWidth = 1.5;  p.Color = clr(8,:);  p.LineStyle = str{1};  p.MarkerSize = 8;  hold on;
Lgd2{4} = sprintf('Gen-D');

hold off;
legend(Lgd2);
grid on;
axis([-Inf, Inf, 0.2, 1.5]);
xlabel('Global Iteration','Interpreter','latex');
ylabel('Training Loss','Interpreter','latex');
set(gca,'linewidth',1.5);
set(gca,'fontsize',15);
set(gca,'Box','off');