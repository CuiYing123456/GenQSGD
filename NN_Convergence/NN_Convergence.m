
clear; close all; clc
% parpool;
%% ============ Part 1: Loading Data =============
pathname_data = '/home/lyc/LYC/2021_Journal/Revision1/NN_Convergence/Data/';

fprintf('Loading Data ...\n')

% load training data set
fid = fopen(strcat(pathname_data, 'train-images.idx3-ubyte'), 'rb');
A_train_img = fread(fid, 4, 'uint32', 'ieee-be');      %锟斤拷取锟侥硷拷头
img_train_num = A_train_img(2);       % 训锟斤拷锟斤拷锟捷硷拷锟斤拷锟斤拷锟斤拷
row_num = A_train_img(3);     %锟斤拷锟斤拷锟斤拷锟斤拷
col_num = A_train_img(4);     %锟斤拷锟斤拷锟斤拷锟斤拷
img_train = fread(fid, [row_num*col_num, img_train_num], 'unsigned char', 'ieee-be');  %锟斤拷取图片
img_train = (img_train/255);        % 锟斤拷一锟斤拷

% data transportation for display and reshape format
img_train = DataTrans( img_train, row_num, col_num );

% load training label set
fid = fopen(strcat(pathname_data, 'train-labels.idx1-ubyte'), 'rb');
A_train_label = fread(fid, 2, 'uint32', 'ieee-be');
lab_train_num = A_train_label(2);       % 训锟斤拷锟斤拷签锟斤拷
lab_train = fread(fid, lab_train_num, 'unsigned char',  'ieee-be');     % 训锟斤拷锟斤拷签

% load test data set
fid = fopen(strcat(pathname_data, 't10k-images.idx3-ubyte'), 'rb');
A_test_img = fread(fid, 4, 'uint32', 'ieee-be');      %锟斤拷取锟侥硷拷头
img_test_num = A_test_img(2);       %锟斤拷锟斤拷锟斤拷锟捷硷拷锟斤拷锟斤拷锟斤拷
img_test = fread(fid, [row_num*col_num, img_test_num], 'unsigned char',  'ieee-be');  %锟斤拷取图片
img_test = (img_test/255);        %锟斤拷一锟斤拷

% data transportation for display and reshape format
img_test = DataTrans( img_test, row_num, col_num );

% load test label set
fid = fopen(strcat(pathname_data, 't10k-labels.idx1-ubyte'), 'rb');
A_test_label = fread(fid, 2, 'uint32', 'ieee-be');
lab_test_num = A_test_label(2);         % 锟斤拷锟皆憋拷签锟斤拷
lab_test = fread(fid, lab_test_num, 'unsigned char',  'ieee-be');
lab_test = lab_test';
%% %%%%%%%%%%%%%% Part 2: Network and Simulation Setting %%%%%%%%%%%%%%%
% network hyperparameter
input_layer_size  = 784;  % 28 x 28 Input Images of Digits
hidden_layer_size = 128;   % number of hidden units
num_labels = 10;          % 10 labels, from 1 to 10
lambda = 0;     % regularization
delta = 0.4;        %convergence threshold

% system hyperparameter
N = 10;     % N workers
local_batch_size = floor(img_train_num / N);     % local batch size -- 锟斤拷锟斤拷60000锟斤拷训锟斤拷锟斤拷锟斤拷平锟斤拷锟街诧拷锟斤拷N锟斤拷worker
batch_size_th = 1;        % local拥锟斤拷锟斤拷小sample size

% simulation hyperparameter
G_th_list = 0.2:0.01:0.25;

% local dataset
local_batch_size = img_train_num/N;
local_img = zeros(row_num * col_num, local_batch_size, N);
local_lab = zeros(local_batch_size, N);
for i = 1:N
    mini_idx = [local_batch_size * (i - 1) + 1: local_batch_size * i];
    local_img(:, :, i) = img_train(:, mini_idx);
    local_lab(:, i) = lab_train(mini_idx);
end

% initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
% initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);
% J_epoch_init = nnCostFunction(initial_Theta1, initial_Theta2, num_labels, img_train, lab_train, lambda);
%
% pred = Predict(initial_Theta1, initial_Theta2, img_train);
% pred = pred';
% Accuracy_epoch_init = mean(double(pred == lab_train)) * 100;
% save(['initial_Theta1.mat'], 'initial_Theta1');
% save(['initial_Theta2.mat'], 'initial_Theta2');
% save(['J_epoch_init.mat'], 'J_epoch_init');
% save(['Accuracy_epoch_init.mat'], 'Accuracy_epoch_init');

load (strcat(pathname_data, 'initial_Theta1.mat'));
load (strcat(pathname_data, 'initial_Theta2.mat'));
load (strcat(pathname_data, 'J_epoch_init.mat'));
load (strcat(pathname_data, 'Accuracy_epoch_init.mat'));

%% ==============  Part 3-1: Simulation for Opt Step Size ==================
s0t = 2^32-1;
s0 = 2^32-1;
snt = 2^32-1;
sn = 2^32-1;
range = 0.3;

% bit_server_worker = {...
%     [32, 32]; [32, 4]; [32, 1]; ...
%     [8, 16]; [8, 4]; [8, 1]; ...
%     [2, 16]; [2, 4]; [2, 1]};
% bit_server_worker = {...
%     [32, 32]; [32, 16]; [32, 8]; [16,32]; [16, 16]; [16, 8]; [8, 32]; [8, 16]; [8, 8]};
% bit_server_worker = {...
%     [8, 8]; [8, 4]; [8, 2]; [4,8]; [4, 4]; [4, 2]; [2, 8]; [2, 4]; [2, 2]};
bit_server_worker = {...
    [4, 4]; [4, 2]; [4, 1]; [2, 4]; [2, 2]; [2, 1]; [1, 4]; [1, 2]; [1, 1]};

% B_list = [8,8,8,8,7,7];
% K_1_list = [1,1,1,1,1,1];
% K_2_list = [1,1,1,1,1,1];
% K_0_list = [1213,1117,1034,960,935,872];      % !!!!!!!300
% lr_list = [0.0515,0.0524,0.0532,0.0541,0.0549,0.0557];
B = 1;
K_1 = 1;
K_2 = 1;
K_0 = 1000;
lr = 0.05;
L = 2;

for i = 1 : length(bit_server_worker)
    bit = bit_server_worker{i};
    sn = 2^bit(2)-1;
    s0 = 2^bit(1)-1;
    %     B = B_list(i);
    %     K_1 = K_1_list(i);
    %     K_2 = K_2_list(i);
    %     K_0 = K_0_list(i);
    %     lr = lr_list(i);
    
    J_step = zeros(K_0, L);
    accur = zeros(K_0, L);
    
    for iter=1:L
        fprintf('\nTraining Neural Network... \n')
        Theta1_global=initial_Theta1;
        Theta2_global=initial_Theta2;
        Theta1_local = zeros(size(Theta1_global));
        Theta2_local = zeros(size(Theta2_global));
        QuantizedGradient1 = zeros([size(Theta1_global), N]);
        QuantizedGradient2 = zeros([size(Theta2_global), N]);
        for k_0=1:K_0
            % model distribution
            for n = 1 : N
                Theta1_local(:, :, n) = Theta1_global;
                Theta2_local(:, :, n) = Theta2_global;
            end
            % sampling
            rand_idx = randperm(img_train_num/N);
            rand_idx = rand_idx(1:B);
            for n = 1 : 5
                Theta1_local = Theta1_global;
                Theta2_local = Theta2_global;
                for k_1 = 1 : K_1
                    X = local_img(:, rand_idx, n);
                    y = local_lab(rand_idx, n);
                    % SGD update
                    [local_Grad1, local_Grad2] = nnGradient(Theta1_local,Theta2_local, input_layer_size, hidden_layer_size, num_labels, X, y, lambda);  %�ݶ�
                    Theta1_local = Theta1_local - lr * local_Grad1;
                    Theta2_local = Theta2_local - lr * local_Grad2;
                end
                [QuantizedGradient1(:, :, n), QuantizedGradient2(:, :, n)] = ...
                    Quantization(Theta1_local - Theta1_global, Theta2_local - Theta2_global, sn);
            end
            for n = 6 : 10
                for k_2 = 1 : K_2
                    X = local_img(:, rand_idx, n);
                    y = local_lab(rand_idx, n);
                    % SGD update
                    [local_Grad1, local_Grad2] = nnGradient(Theta1_local,Theta2_local, input_layer_size, hidden_layer_size, num_labels, X, y, lambda);  %�ݶ�
                    Theta1_local = Theta1_local - lr * local_Grad1;
                    Theta2_local = Theta2_local - lr * local_Grad2;
                end
                [QuantizedGradient1(:, :, n), QuantizedGradient2(:, :, n)] = ...
                    Quantization(Theta1_local - Theta1_global, Theta2_local - Theta2_global, sn);
            end
            % model update
            AverageGradient1 = mean(QuantizedGradient1, 3);
            AverageGradient2 = mean(QuantizedGradient2, 3);
            [QuantizedAverage1, QuantizedAverage2] = Quantization(AverageGradient1, AverageGradient2, s0);
            Theta1_global = Theta1_global + QuantizedAverage1;
            Theta2_global = Theta2_global + QuantizedAverage2;
            
            J = nnCostFunction(Theta1_global, Theta2_global, num_labels, img_train, lab_train, lambda);
            J_step(k_0,iter)=J;
            pred = Predict(Theta1_global, Theta2_global, img_test);
            accur(k_0,iter)=mean(double(pred == lab_test)) * 100;
            fprintf('\nTest accuracy of step %d of iteration %d for i=%d (SGD): %f', k_0, iter, i,...
                mean(double(pred == lab_test)) * 100);
        end
        
        a_mean=mean(accur,2);
        c_mean=mean(J_step,2);
        save(sprintf('./Data/acc_b0%dbn%d.mat', bit(1), bit(2)), 'a_mean');
        save(sprintf('./Data/cost_b0%dbn%d.mat', bit(1), bit(2)), 'c_mean');
    end
end

% %% =============== Part 4: Plot ==========
% clr=[[1 0 0];[1 0 1];[0 0 1];[0 0.5 0.5];[0.2 0.5 0.3];[0 0 0]];
% str=[' -';' :';'--'; '-.'];
% mkr = ['o';'+';'s';'^';'*';'d'];
% 
% figure(1);
% yyaxis left;        % N-left-time
% Lgd{1} = 'accuracy (A)';
% load(strcat(pathname_data, 'Acc_rcd_opt.mat'));
% Acc_max = max(Acc_rcd_opt);
% Acc_min = min(Acc_rcd_opt);
% p1 = plot(G_th_list, Acc_rcd_opt);
% p1.Color = 'r';
% p1.LineStyle = str(1,:);
% p1.Marker = mkr(1,:);
% p1.MarkerSize = 10;
% p1.LineWidth = 1.5;
% hold on;
% Lgd{2} = 'accuracy (C)';
% load(strcat(pathname_data, 'Acc_rcd_cons.mat'));
% Acc_max = max( [Acc_max, Acc_rcd_cons] );
% Acc_min = min( [Acc_min, Acc_rcd_cons] );
% p1 = plot(G_th_list, Acc_rcd_cons);
% p1.Color = 'r';
% p1.LineStyle = str(1,:);
% p1.Marker = mkr(2,:);
% p1.MarkerSize = 10;
% p1.LineWidth = 1.5;
% hold on;
% Lgd{3} = 'accuracy (E)';
% load(strcat(pathname_data, 'Acc_rcd_exp.mat'));
% Acc_max = max( [Acc_max, Acc_rcd_exp] );
% Acc_min = min( [Acc_min, Acc_rcd_exp] );
% p1 = plot(G_th_list, Acc_rcd_exp);
% p1.Color = 'r';
% p1.LineStyle = str(1,:);
% p1.Marker = mkr(3,:);
% p1.MarkerSize = 10;
% p1.LineWidth = 1.5;
% hold on;
% Lgd{4} = 'accuracy (D)';
% load(strcat(pathname_data, 'Acc_rcd_dim.mat'));
% Acc_max = max( [Acc_max, Acc_rcd_dim] );
% Acc_min = min( [Acc_min, Acc_rcd_dim] );
% p1 = plot(G_th_list, Acc_rcd_dim);
% p1.Color = 'r';
% p1.LineStyle = str(1,:);
% p1.Marker = mkr(4,:);
% p1.MarkerSize = 10;
% p1.LineWidth = 1.5;
% hold off;
% 
% axis([min(G_th_list), max(G_th_list), 0.95*Acc_min 1.02*Acc_max]);
% ylabel('Test Accuracy ($\%$)','interpreter','latex','fontsize',12);
% set(gca,'ycolor','r');
% set(gca,'linewidth',1.5);
% set(gca,'FontSize',12);
% 
% 
% 
% yyaxis right;        % N-left-time
% Lgd{5} = 'loss (A)';
% load(strcat(pathname_data, 'Loss_rcd_opt.mat'));
% Loss_max = max(Loss_rcd_opt);
% Loss_min = min(Loss_rcd_opt);
% p1 = plot(G_th_list, Loss_rcd_opt);
% p1.Color = 'b';
% p1.LineStyle = str(2,:);
% p1.Marker = mkr(1,:);
% p1.MarkerSize = 10;
% p1.LineWidth = 1.5;
% hold on;
% Lgd{6} = 'loss (C)';
% load(strcat(pathname_data, 'Loss_rcd_cons.mat'));
% Loss_max = max( [Loss_rcd_cons, Loss_max] );
% Loss_min = min( [Loss_rcd_cons, Loss_min] );
% p1 = plot(G_th_list, Loss_rcd_cons);
% p1.Color = 'b';
% p1.LineStyle = str(2,:);
% p1.Marker = mkr(2,:);
% p1.MarkerSize = 10;
% p1.LineWidth = 1.5;
% hold on;
% Lgd{7} = 'loss (E)';
% load(strcat(pathname_data, 'Loss_rcd_exp.mat'));
% Loss_max = max( [Loss_rcd_exp, Loss_max] );
% Loss_min = min( [Loss_rcd_exp, Loss_min] );
% p1 = plot(G_th_list, Loss_rcd_exp);
% p1.Color = 'b';
% p1.LineStyle = str(2,:);
% p1.Marker = mkr(3,:);
% p1.MarkerSize = 10;
% p1.LineWidth = 1.5;
% hold on;
% Lgd{8} = 'loss (D)';
% load(strcat(pathname_data, 'Loss_rcd_dim.mat'));
% Loss_max = max( [Loss_rcd_dim, Loss_max] );
% Loss_min = min( [Loss_rcd_dim, Loss_min] );
% p1 = plot(G_th_list, Loss_rcd_dim);
% p1.Color = 'b';
% p1.LineStyle = str(2,:);
% p1.Marker = mkr(4,:);
% p1.MarkerSize = 10;
% p1.LineWidth = 1.5;
% hold off;
% 
% axis([ min(G_th_list), max(G_th_list), 0.9*Loss_min 1.5*Loss_max]);
% ax = gca;
% ax.Box = 'off';
% ylabel('Training Loss','interpreter','latex','fontsize',12);
% xlabel('Convergence Error Limit $C_{\max}$','interpreter','latex','fontsize',12);
% set(gca,'ycolor','b');
% set(gca,'linewidth',1.2);
% set(gca,'FontSize',12);
% columnlegend(2, Lgd);
% 
% grid on

		

