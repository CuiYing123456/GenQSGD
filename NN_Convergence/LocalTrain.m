function [ Theta1, Theta2] = LocalTrain( Theta1_avg, Theta2_avg, ...
                                                                 input_layer_size, hidden_layer_size, num_labels, ...
                                                                 local_img, local_lab, B, step, lr, lambda, k)

Theta1 = Theta1_avg;
Theta2 = Theta2_avg;
% [Theta1_grad_init, Theta2_grad_init] = nnGradient(Theta1,Theta2, input_layer_size, hidden_layer_size, num_labels, local_img, local_lab, lambda);
batch_size = size(local_img,2);
local_epoch_size = step * B;
used_size = (k-1) * local_epoch_size;       %local����sample��
start_idx = mod(used_size, batch_size) + 1;        % ����ѵ����ʼ������
reuse_batch_num = floor(local_epoch_size / batch_size);     % ����ѵ����Ҫ����batch�Ĵ���
sel_idx = [];

for t = 1 : step

    sel_idx = randi( batch_size, 1, B);
    X = local_img(:, sel_idx);      % stochastic minibatch sample set
    y = local_lab(sel_idx);     % stochastic minibatch label set
    [Theta1_grad, Theta2_grad] = nnGradient(Theta1,Theta2, input_layer_size, hidden_layer_size, num_labels, X, y, lambda);  %�ݶ�
    Theta1 = Theta1 - lr * Theta1_grad;     %gradient descent
    Theta2 = Theta2 - lr * Theta2_grad;
    
end




