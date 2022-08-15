function [J_opt, J_epoch, Acc_opt, Accuracy_epoch, epoch_end, IsOverflow] ...
    = GlobalTrain_dim( num_labels, img_train, lab_train, img_test, lab_test, ...
    img_train_num, row_num, col_num, local_batch_size, ...
    input_layer_size, hidden_layer_size, ...
    E, N, step1, step2, B, gamma_D, rho_D, lambda,...
    initial_Theta1,initial_Theta2 , J_epoch_init , Accuracy_epoch_init)
fprintf('\nGlobal training B = %d... \n', B)

% initialize data set
global_rand = randperm(img_train_num);
img_train_rand = img_train(:, global_rand);        %ѵ����������
lab_train_rand = lab_train(global_rand);        %ѵ����ǩ����
local_img = zeros(row_num * col_num, local_batch_size, N);      % local������
local_lab = zeros(local_batch_size, N);       %local��ǩ��

for i = 1:N
    mini_idx = [local_batch_size * (i - 1) + 1: local_batch_size * i];
    local_img(:, :, i) = img_train_rand(:, mini_idx);
    local_lab(:, i) = lab_train_rand(mini_idx);
end
fprintf('\nData allocation finished!\n');
Theta1_avg = initial_Theta1;
Theta2_avg = initial_Theta2;

% Initialize performance
%J_epoch = nnCostFunction(Theta1_avg, Theta2_avg, num_labels, img_train, lab_train, lambda);

J_epoch = J_epoch_init;
Accuracy_epoch = Accuracy_epoch_init;
IsOverflow = 0;
J_opt = J_epoch_init;


% NN training
for k=1 : E
    fprintf('\n==== epoch %d start ====\n',k)
    
    Theta1_accu = zeros(size(initial_Theta1));
    Theta2_accu = zeros(size(initial_Theta2));
    lr = rho_D*gamma_D/(k+rho_D);
    % local training
    for i=1 : floor(N/2)
        [local_Theta1, local_Theta2] = LocalTrain( Theta1_avg, Theta2_avg, ...
            input_layer_size, hidden_layer_size, num_labels, ...
            local_img(:, :, i), local_lab(:, i), B, step1, lr, lambda, k);
        Theta1_accu = Theta1_accu + local_Theta1;
        Theta2_accu = Theta2_accu + local_Theta2;
    end
    for i= floor(N/2) + 1 : N
        [local_Theta1, local_Theta2] = LocalTrain( Theta1_avg, Theta2_avg, ...
            input_layer_size, hidden_layer_size, num_labels, ...
            local_img(:, :, i), local_lab(:, i), B, step2, lr, lambda, k);
        Theta1_accu = Theta1_accu + local_Theta1;
        Theta2_accu = Theta2_accu + local_Theta2;
    end
    % global averaging
    Theta1_avg = Theta1_accu / N;

    % global cost function
    J_tmp = nnCostFunction(Theta1_avg, Theta2_avg, num_labels, img_train, lab_train, lambda);
    if isnan(J_tmp)
        fprintf('\nGlobal training overflow, last cost  = %f, last accuracy = %d\n', J_epoch(k), Accuracy_epoch(k));
        IsOverflow = 1;
        epoch_end = k;
        return
    end
    J_epoch = [J_epoch, J_tmp];
    fprintf('Loss: %f \n',J_tmp );
    
    % test accuracy
    pred = Predict(Theta1_avg, Theta2_avg, img_test);
    pred = pred';
    A_tmp = mean(double(pred == lab_test)) * 100;
    Accuracy_epoch = [Accuracy_epoch, A_tmp];
    fprintf('\nTraining Set Accuracy: %2.2f%%\n', A_tmp);    %׼ȷ��

    if k==E %���ܴﲻ��������ֵ
        epoch_end = k;
        J_opt = J_epoch(k+1);
        Acc_opt = Accuracy_epoch(k+1);
        fprintf('\nlast cost  = %f, last accuracy = %d, epoch_end = %d\n', ...
            J_epoch(k + 1), Accuracy_epoch(k + 1), epoch_end);
        return
    end
end

end