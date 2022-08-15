clr = [[1, 0, 0]; [1, 0.5, 0]; [1, 0.8, 0.2]; 
    [0, 0.5, 0]; [0.2, 0.8, 0]; [0.6, 1, 0];
    [0, 0.2, 1]; [0, 0.6, 1]; [0, 1, 1]];
str = {'-'; '--'; '-.'; ':'};

%% =========== Plot Test Accuracy ===========
Lgd1 = {};
figure(1);
% bit_server_worker = {...
%     [32, 32]; [32, 4]; [32, 1]; ...
%     [8, 16]; [8, 4]; [8, 1]; ...
%     [2, 16]; [2, 4]; [2, 1]};
% bit_server_worker = {...
%     [32, 32]; [32, 16]; [32, 8]; [16,32]; [16, 16]; [16, 8]; [8, 32]; [8, 16]; [8, 8]};
bit_server_worker = {...
    [8, 8]; [8, 4]; [8, 2]; [4,8]; [4, 4]; [4, 2]; [2, 8]; [2, 4]; [2, 2]};
% bit_server_worker = {...
%     [4, 4]; [4, 2]; [4, 1]; [2, 4]; [2, 2]; [2, 1]; [1, 4]; [1, 2]; [1, 1]};

bit = bit_server_worker{1};
load(sprintf('./Data/acc_b0%dbn%d.mat', bit(1), bit(2)), 'a_mean');
p = plot(a_mean);  p.LineWidth = 1.2;  p.Color = clr(1,:);  p.LineStyle = str{1};  hold on;
Lgd1{1} = sprintf('b0=%d,bn=%d', bit(1), bit(2));
bit = bit_server_worker{2};
load(sprintf('./Data/acc_b0%dbn%d.mat', bit(1), bit(2)), 'a_mean');
p = plot(a_mean);  p.LineWidth = 1.2;  p.Color = clr(2,:);  p.LineStyle = str{1};  hold on;
Lgd1{2} = sprintf('b0=%d,bn=%d', bit(1), bit(2));
bit = bit_server_worker{3};
load(sprintf('./Data/acc_b0%dbn%d.mat', bit(1), bit(2)), 'a_mean');
p = plot(a_mean);  p.LineWidth = 1.2;  p.Color = clr(3,:);  p.LineStyle = str{1};  hold on;
Lgd1{3} = sprintf('b0=%d,bn=%d', bit(1), bit(2));

bit = bit_server_worker{4};
load(sprintf('./Data/acc_b0%dbn%d.mat', bit(1), bit(2)), 'a_mean');
p = plot(a_mean);  p.LineWidth = 1.2;  p.Color = clr(4,:);  p.LineStyle = str{1};  hold on;
Lgd1{4} = sprintf('b0=%d,bn=%d', bit(1), bit(2));
bit = bit_server_worker{5};
load(sprintf('./Data/acc_b0%dbn%d.mat', bit(1), bit(2)), 'a_mean');
p = plot(a_mean);  p.LineWidth = 1.2;  p.Color = clr(5,:);  p.LineStyle = str{1};  hold on;
Lgd1{5} = sprintf('b0=%d,bn=%d', bit(1), bit(2));
bit = bit_server_worker{6};
load(sprintf('./Data/acc_b0%dbn%d.mat', bit(1), bit(2)), 'a_mean');
p = plot(a_mean);  p.LineWidth = 1.2;  p.Color = clr(6,:);  p.LineStyle = str{1};  hold on;
Lgd1{6} = sprintf('b0=%d,bn=%d', bit(1), bit(2));

bit = bit_server_worker{7};
load(sprintf('./Data/acc_b0%dbn%d.mat', bit(1), bit(2)), 'a_mean');
p = plot(a_mean);  p.LineWidth = 1.2;  p.Color = clr(7,:);  p.LineStyle = str{1};  hold on;
Lgd1{7} = sprintf('b0=%d,bn=%d', bit(1), bit(2));
bit = bit_server_worker{8};
load(sprintf('./Data/acc_b0%dbn%d.mat', bit(1), bit(2)), 'a_mean');
p = plot(a_mean);  p.LineWidth = 1.2;  p.Color = clr(8,:);  p.LineStyle = str{1};  hold on;
Lgd1{8} = sprintf('b0=%d,bn=%d', bit(1), bit(2));
bit = bit_server_worker{9};
load(sprintf('./Data/acc_b0%dbn%d.mat', bit(1), bit(2)), 'a_mean');
p = plot(a_mean);  p.LineWidth = 1.2;  p.Color = clr(9,:);  p.LineStyle = str{1};  hold on;
Lgd1{9} = sprintf('b0=%d,bn=%d', bit(1), bit(2));


hold off;
legend(Lgd1);
grid on;
axis([-Inf, Inf, 50, 100]);
xlabel('Iteration');
ylabel('Test accuracy');

%% =========== Plot Cost Function ===========
Lgd2 = {};
figure(2);
cnt = 1;

load(sprintf('./SaveData/SCc1_B10.mat'));
p = plot(c_mean1);  p.LineWidth = 1.2;  p.Color = clr(1,:);  p.LineStyle = str{1};  hold on;
Lgd2{1} = sprintf('Alg.1,B=10');
load(sprintf('./SaveData/SCc1_B100.mat'));
p = plot(c_mean1);  p.LineWidth = 1.2;  p.Color = clr(2,:);  p.LineStyle = str{1};  hold on;
Lgd2{2} = sprintf('Alg.1,B=100');
load(sprintf('./SaveData/SCc1_B6000.mat'));
p = plot(c_mean1);  p.LineWidth = 1.2;  p.Color = clr(3,:);  p.LineStyle = str{1};  hold on;
Lgd2{3} = sprintf('Alg.1,B=6000');

% for i = 1 : length(B)
%     load(sprintf('./SaveData/SCc1_B%d', B(i)));
%     p = plot(c_mean1, 'LineWidth', 1.5);  p.Color = clr(i,:);  p.LineStyle = str{1};  hold on;
% end
% 
% for i = 1 : length(B)
%     Lgd2{cnt} = strcat('Alg.1, B=', num2str(B(i)));
%     cnt = cnt + 1;
% end


hold off;
legend(Lgd2);
grid on;
axis([-Inf, Inf, 0, 1.5]);
xlabel('Iteration');
ylabel('Training cost');