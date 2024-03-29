function [ energy ] = E_Q_FedAvg( N, s_n, s_0, D, alpha, ...
    C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
    Q, gamma, G, L, sigma, C_max, T_max )

q_0 = min(D/s_0^2, sqrt(D)/s_0);
q_n = min([D./(s_n.^2), sqrt(D)./s_n], [], 2);
M_n = 32 + D + D*log2(s_n);
M_0 = 32 + D + D*log2(s_0);

c_1 = 2*N*Q;
c_2 = 4*G^2*L^2;
c_3 = L*sigma^2/N;
c_4 = 2*L*G^2;
q_s0sn = q_0+q_n+q_0*q_n;

%% =============== Solution: FedAvg =================
% ------------- Initialize --------------
K_t = 0;
K_0_t = 0;
B_t = 0;

K_init_list = [1:50, 60:20:1000, 1500:500:10000];
K_0_init_list = [1:50, 100, 1000, 10000];
minibatch_sz = 100;
batch_num_list = [1:2:100];
% % ---------- Check -------------
is_feasible = 0;
for i1 = 1 : length(K_init_list)
    K = K_init_list(i1);
    for i2 = 1 : length(K_0_init_list)
        K_0 = K_0_init_list(i2);
        for i3 = 1 : length(batch_num_list)
            batch_num = batch_num_list(i3);
            B = batch_num*minibatch_sz/K;
            cons3 = isempty(   find( ((B*K*max(C_n./F_n)+C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 <= T_max)==0, 1 )   );
            cons4 = isempty(   find( (c_1/(gamma*K_0*N*K) + c_2*gamma^2*K^2 ...
                + c_3*gamma/B + c_4*gamma*sum(q_s0sn)*K/N <= C_max)==0, 1)   );

            if cons3 && cons4
                is_feasible = 1;    fprintf('Feasible initial point found!\n');   break;
            end
        end
        if is_feasible
            break;
        end
    end
    if is_feasible
        break;
    end
end

% ------------- Algorithm --------------
energy_rcd = [];
norm_rcd = [];
while norm([K_t; K_0_t; B_t]-[K; K_0; B], 2)>0.1
    K_t = K;
    K_0_t = K_0;
    B_t = B;

    cvx_begin gp quiet
    cvx_solver SeDuMi
    variables K K_0 B
    minimize (   K_0 * (B*K*sum(alpha*C_n.*F_n.^2) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0)  )
    subject to
    (C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 + B*K_0*K*max(C_n./F_n) <= T_max;
    c_1/(gamma*K_0*N*K) + c_2*gamma^2*K^2 ...
        + c_3*gamma/B + c_4*gamma*sum(q_s0sn)*K/N <= C_max;
    K * B <= batch_num*minibatch_sz;
    K * B >= batch_num*minibatch_sz;
    cvx_end

    energy_rcd = [energy_rcd, K_0 * (B*K*sum(alpha*C_n.*F_n.^2) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0)];
    norm_rcd = [norm_rcd, norm([K_t; K_0_t; B_t]-[K; K_0; B], 2)];
end
energy = K_0 * (B*K*sum(alpha*C_n.*F_n.^2) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);

K(K<1) = 1;
B(B<1) = 1;
K_0(K_0<1) = 1;

K_tmp = [ ceil( K ), floor( K ) ]; B_tmp = [ ceil(B), floor(B) ]; K_0_tmp = [ ceil(K_0), floor(K_0) ];   rcv_error = 10000000000;

for i = 1 : size(K_tmp, 2)
    K = K_tmp(:, i);
    for j = 1 : length(K_0_tmp)
        K_0 = K_0_tmp(j);
        for k = 1 : length(B_tmp)
            B = B_tmp(k);
            energy_tmp = K_0 * (B*K*sum(alpha*C_n.*F_n.^2) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);
            if abs(energy_tmp-energy_rcd(end)) < rcv_error
                energy = energy_tmp;
                rcv_error = abs(energy_tmp-energy_rcd(end));
            end
        end
    end
end
