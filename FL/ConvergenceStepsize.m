%% Start
clear; close all; clc
%% ============ Part 1: Loading Data =============
pathname_data = './Data/';

fprintf('Loading Data ...\n')

% load training data set
fid = fopen(strcat(pathname_data, 'train-images.idx3-ubyte'), 'rb');
A_train_img = fread(fid, 4, 'uint32', 'ieee-be');     
img_train_num = A_train_img(2);     
row_num = A_train_img(3);   
col_num = A_train_img(4);   
img_train = fread(fid, [row_num*col_num, img_train_num], 'unsigned char', 'ieee-be'); 
img_train = (img_train/255);     

% data transportation for display and reshape format
img_train = DataTrans( img_train, row_num, col_num );

% load training label set
fid = fopen(strcat(pathname_data, 'train-labels.idx1-ubyte'), 'rb');
A_train_label = fread(fid, 2, 'uint32', 'ieee-be');
lab_train_num = A_train_label(2);      
lab_train = fread(fid, lab_train_num, 'unsigned char',  'ieee-be');     

% load test data set
fid = fopen(strcat(pathname_data, 't10k-images.idx3-ubyte'), 'rb');
A_test_img = fread(fid, 4, 'uint32', 'ieee-be');    
img_test_num = A_test_img(2);   
img_test = fread(fid, [row_num*col_num, img_test_num], 'unsigned char',  'ieee-be'); 
img_test = (img_test/255);       

% data transportation for display and reshape format
img_test = DataTrans( img_test, row_num, col_num );

% load test label set
fid = fopen(strcat(pathname_data, 't10k-labels.idx1-ubyte'), 'rb');
A_test_label = fread(fid, 2, 'uint32', 'ieee-be');
lab_test_num = A_test_label(2);       
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
batch_size_th = 1;        

sn_mean = 2^14;
sn_ratio = 10;
sn_min = 2*sn_mean/(1+sn_ratio);
sn_max = 2*sn_ratio*sn_mean/(1+sn_ratio);
sn1 = sn_max;
sn2 = sn_min;
s0 = 2^32;

% local dataset
local_batch_size = img_train_num/N;
local_img = zeros(row_num * col_num, local_batch_size, N);
local_lab = zeros(local_batch_size, N);
for i = 1:N
    mini_idx = [local_batch_size * (i - 1) + 1: local_batch_size * i];
    local_img(:, :, i) = img_train(:, mini_idx);
    local_lab(:, i) = lab_train(mini_idx);
end
% simulation hyperparameter
G_th_list = 0.2:0.01:0.25;
L = 10;

%% ==============  Part 3-1: Simulation for Opt Step Size ==================
B_list = [8,8,8,7,7,7];
K_1_list = [1,1,1,1,1,1];
K_2_list = [1,1,1,1,1,1];
K_0_list = [1213,1117,1034,1004,935,872];   
lr_list = [0.0515,0.0519,0.0528,0.0537,0.0545,0.0554];

for i = 1 : 1
    B = B_list(i);  K_1 = K_1_list(i);  K_2 = K_2_list(i);  K_0 = K_0_list(i);  lr = lr_list(i);
    J_step = zeros(K_0, L);  accur = zeros(K_0, L);
    for iter=1:L
        % initialization
        initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
        initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);
        % training
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
                    [local_Grad1, local_Grad2] = nnGradient(Theta1_local,Theta2_local, input_layer_size, hidden_layer_size, num_labels, X, y, lambda);                     Theta1_local = Theta1_local - lr * local_Grad1;
                    Theta2_local = Theta2_local - lr * local_Grad2;
                end
                [QuantizedGradient1(:, :, n), QuantizedGradient2(:, :, n)] = ...
                    Quantization(Theta1_local - Theta1_global, Theta2_local - Theta2_global, sn1);
            end
            for n = 6 : 10
                for k_2 = 1 : K_2
                    X = local_img(:, rand_idx, n);
                    y = local_lab(rand_idx, n);
                    % SGD update
                    [local_Grad1, local_Grad2] = nnGradient(Theta1_local,Theta2_local, input_layer_size, hidden_layer_size, num_labels, X, y, lambda);  
                    Theta1_local = Theta1_local - lr * local_Grad1;
                    Theta2_local = Theta2_local - lr * local_Grad2;
                end
                [QuantizedGradient1(:, :, n), QuantizedGradient2(:, :, n)] = ...
                    Quantization(Theta1_local - Theta1_global, Theta2_local - Theta2_global, sn2);
            end
            % model update
            AverageGradient1 = mean(QuantizedGradient1, 3);
            AverageGradient2 = mean(QuantizedGradient2, 3);
            [QuantizedAverage1, QuantizedAverage2] = Quantization(AverageGradient1, AverageGradient2, s0);
            Theta1_global = Theta1_global + QuantizedAverage1;
            Theta2_global = Theta2_global + QuantizedAverage2;
            
            J_step(k_0,iter) = nnCostFunction(Theta1_global, Theta2_global, num_labels, img_train, lab_train, lambda);
            pred = Predict(Theta1_global, Theta2_global, img_test);
            accur(k_0,iter)=mean(double(pred == lab_test)) * 100;
            fprintf('\nTest accuracy of step %d of iteration %d (opt): %f', k_0, iter,...
                mean(double(pred == lab_test)) * 100);
        end
        
        a_mean=mean(accur,2);
        c_mean=mean(J_step,2);
        save(sprintf('./Data/acc_opt.mat'), 'a_mean');
        save(sprintf('./Data/cost_opt.mat'), 'c_mean');
    end  
end

%% ==============  Part 3-2: Simulation for Constant Step Size ==================
gamma_C = 0.01;
B_list = [2,2,2,2,3,3];
K_1_list = [3,3,3,3,3,3];
K_2_list = [3,3,3,3,3,3];
K_0_list = [1357,1255,1168,1092,935,886];

for i = 1 : 1
    B = B_list(i);  K_1 = K_1_list(i);  K_2 = K_2_list(i);  K_0 = K_0_list(i);  lr = gamma_C;
    J_step = zeros(K_0, L);  accur = zeros(K_0, L);
    for iter=1:L
        % initialization
        initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
        initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);
        % training
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
                    [local_Grad1, local_Grad2] = nnGradient(Theta1_local,Theta2_local, input_layer_size, hidden_layer_size, num_labels, X, y, lambda); 
                    Theta1_local = Theta1_local - lr * local_Grad1;
                    Theta2_local = Theta2_local - lr * local_Grad2;
                end
                [QuantizedGradient1(:, :, n), QuantizedGradient2(:, :, n)] = ...
                    Quantization(Theta1_local - Theta1_global, Theta2_local - Theta2_global, sn1);
            end
            for n = 6 : 10
                for k_2 = 1 : K_2
                    X = local_img(:, rand_idx, n);
                    y = local_lab(rand_idx, n);
                    % SGD update
                    [local_Grad1, local_Grad2] = nnGradient(Theta1_local,Theta2_local, input_layer_size, hidden_layer_size, num_labels, X, y, lambda); 
                    Theta1_local = Theta1_local - lr * local_Grad1;
                    Theta2_local = Theta2_local - lr * local_Grad2;
                end
                [QuantizedGradient1(:, :, n), QuantizedGradient2(:, :, n)] = ...
                    Quantization(Theta1_local - Theta1_global, Theta2_local - Theta2_global, sn2);
            end
            % model update
            AverageGradient1 = mean(QuantizedGradient1, 3);
            AverageGradient2 = mean(QuantizedGradient2, 3);
            [QuantizedAverage1, QuantizedAverage2] = Quantization(AverageGradient1, AverageGradient2, s0);
            Theta1_global = Theta1_global + QuantizedAverage1;
            Theta2_global = Theta2_global + QuantizedAverage2;
            
            J_step(k_0,iter) = nnCostFunction(Theta1_global, Theta2_global, num_labels, img_train, lab_train, lambda);
            pred = Predict(Theta1_global, Theta2_global, img_test);
            accur(k_0,iter)=mean(double(pred == lab_test)) * 100;
            fprintf('\nTest accuracy of step %d of iteration %d (cons): %f', k_0, iter,...
                mean(double(pred == lab_test)) * 100);
        end
        
        a_mean=mean(accur,2);
        c_mean=mean(J_step,2);
        save(sprintf('./Data/acc_cons.mat'), 'a_mean');
        save(sprintf('./Data/cost_cons.mat'), 'c_mean');
    end 
end

%% ==============  Part 3-3: Simulation for Exponential Step Size ==================
Acc_rcd_exp = [];
Loss_rcd_exp = [];
gamma_E = 0.02;
rho_E = 0.9995;
B_list = [3,3,3,3,3,3];
K_1_list = [2,2,2,2,2,2];
K_2_list = [2,2,2,2,2,2];
K_0_list = [1428,1306,1200,1111,1032,963]; 

for i = 1 : 1
    B = B_list(i);  K_1 = K_1_list(i);  K_2 = K_2_list(i);  K_0 = K_0_list(i);
    J_step = zeros(K_0, L);  accur = zeros(K_0, L);
    for iter=1:L
        % initialization
        initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
        initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);
        % training
        fprintf('\nTraining Neural Network... \n')
        Theta1_global=initial_Theta1;
        Theta2_global=initial_Theta2;
        Theta1_local = zeros(size(Theta1_global));
        Theta2_local = zeros(size(Theta2_global));
        QuantizedGradient1 = zeros([size(Theta1_global), N]);
        QuantizedGradient2 = zeros([size(Theta2_global), N]);
        for k_0=1:K_0
            lr = rho_E^(k_0-1)*gamma_E;
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
                    [local_Grad1, local_Grad2] = nnGradient(Theta1_local,Theta2_local, input_layer_size, hidden_layer_size, num_labels, X, y, lambda); 
                    Theta1_local = Theta1_local - lr * local_Grad1;
                    Theta2_local = Theta2_local - lr * local_Grad2;
                end
                [QuantizedGradient1(:, :, n), QuantizedGradient2(:, :, n)] = ...
                    Quantization(Theta1_local - Theta1_global, Theta2_local - Theta2_global, sn1);
            end
            for n = 6 : 10
                for k_2 = 1 : K_2
                    X = local_img(:, rand_idx, n);
                    y = local_lab(rand_idx, n);
                    % SGD update
                    [local_Grad1, local_Grad2] = nnGradient(Theta1_local,Theta2_local, input_layer_size, hidden_layer_size, num_labels, X, y, lambda); 
                    Theta1_local = Theta1_local - lr * local_Grad1;
                    Theta2_local = Theta2_local - lr * local_Grad2;
                end
                [QuantizedGradient1(:, :, n), QuantizedGradient2(:, :, n)] = ...
                    Quantization(Theta1_local - Theta1_global, Theta2_local - Theta2_global, sn2);
            end
            % model update
            AverageGradient1 = mean(QuantizedGradient1, 3);
            AverageGradient2 = mean(QuantizedGradient2, 3);
            [QuantizedAverage1, QuantizedAverage2] = Quantization(AverageGradient1, AverageGradient2, s0);
            Theta1_global = Theta1_global + QuantizedAverage1;
            Theta2_global = Theta2_global + QuantizedAverage2;
            
            J_step(k_0,iter) = nnCostFunction(Theta1_global, Theta2_global, num_labels, img_train, lab_train, lambda);
            pred = Predict(Theta1_global, Theta2_global, img_test);
            accur(k_0,iter)=mean(double(pred == lab_test)) * 100;
            fprintf('\nTest accuracy of step %d of iteration %d (exp): %f', k_0, iter,...
                mean(double(pred == lab_test)) * 100);
        end
        
        a_mean=mean(accur,2);
        c_mean=mean(J_step,2);
        save(sprintf('./Data/acc_exp.mat'), 'a_mean');
        save(sprintf('./Data/cost_exp.mat'), 'c_mean');
    end 
end

%% ==============  Part 3-4: Simulation for Diminishing Step Size ==================
Acc_rcd_dim = [];
Loss_rcd_dim = [];
gamma_D = 0.02;
rho_D = 600;
B_list = [3,3,3,3,3,3];
K_1_list = [3,3,3,3,3,3];
K_2_list = [3,3,3,3,3,3];
K_0_list = [1631,1495,1379,1280,1192,1115]; 
for i = 1 : length(B_list)
    B = B_list(i);  K_1 = K_1_list(i);  K_2 = K_2_list(i);  K_0 = K_0_list(i);
    J_step = zeros(K_0, L);  accur = zeros(K_0, L);
    for iter=1:L
        % initialization
        initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
        initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);
        % training
        fprintf('\nTraining Neural Network... \n')
        Theta1_global=initial_Theta1;
        Theta2_global=initial_Theta2;
        Theta1_local = zeros(size(Theta1_global));
        Theta2_local = zeros(size(Theta2_global));
        QuantizedGradient1 = zeros([size(Theta1_global), N]);
        QuantizedGradient2 = zeros([size(Theta2_global), N]);
        for k_0=1:K_0
            lr = rho_D*gamma_D/(k_0+rho_D);
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
                    [local_Grad1, local_Grad2] = nnGradient(Theta1_local,Theta2_local, input_layer_size, hidden_layer_size, num_labels, X, y, lambda); 
                    Theta1_local = Theta1_local - lr * local_Grad1;
                    Theta2_local = Theta2_local - lr * local_Grad2;
                end
                [QuantizedGradient1(:, :, n), QuantizedGradient2(:, :, n)] = ...
                    Quantization(Theta1_local - Theta1_global, Theta2_local - Theta2_global, sn1);
            end
            for n = 6 : 10
                for k_2 = 1 : K_2
                    X = local_img(:, rand_idx, n);
                    y = local_lab(rand_idx, n);
                    % SGD update
                    [local_Grad1, local_Grad2] = nnGradient(Theta1_local,Theta2_local, input_layer_size, hidden_layer_size, num_labels, X, y, lambda);  
                    Theta1_local = Theta1_local - lr * local_Grad1;
                    Theta2_local = Theta2_local - lr * local_Grad2;
                end
                [QuantizedGradient1(:, :, n), QuantizedGradient2(:, :, n)] = ...
                    Quantization(Theta1_local - Theta1_global, Theta2_local - Theta2_global, sn2);
            end
            % model update
            AverageGradient1 = mean(QuantizedGradient1, 3);
            AverageGradient2 = mean(QuantizedGradient2, 3);
            [QuantizedAverage1, QuantizedAverage2] = Quantization(AverageGradient1, AverageGradient2, s0);
            Theta1_global = Theta1_global + QuantizedAverage1;
            Theta2_global = Theta2_global + QuantizedAverage2;
            
            J_step(k_0,iter) = nnCostFunction(Theta1_global, Theta2_global, num_labels, img_train, lab_train, lambda);
            pred = Predict(Theta1_global, Theta2_global, img_test);
            accur(k_0,iter)=mean(double(pred == lab_test)) * 100;
            fprintf('\nTest accuracy of step %d of iteration %d (dim): %f', k_0, iter,...
                mean(double(pred == lab_test)) * 100);
        end
        
        a_mean=mean(accur,2);
        c_mean=mean(J_step,2);
        save(sprintf('./Data/acc_dim.mat'), 'a_mean');
        save(sprintf('./Data/cost_dim.mat'), 'c_mean');
    end
end
