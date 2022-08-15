function [ energy] = E_Q_PR( N, s_n, s_0, D, alpha, ...
    C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
    Q, gamma, G, L, sigma, C_max, T_max )

energy = Inf;

q_0 = min(D/s_0^2, sqrt(D)/s_0);
q_n = min([D./(s_n.^2), sqrt(D)./s_n], [], 2);
M_n = 32 + D + D*log2(s_n);
M_0 = 32 + D + D*log2(s_0);

c_1 = 2*N*Q;
c_2 = 4*G^2*L^2;
c_3 = L*sigma^2/N;
c_4 = 2*L*G^2;
q_s0sn = q_0+q_n+q_0*q_n;

%% =============== Solution: Proposed =================
% ------------- Initialize --------------
K_t = zeros(N, 1);
K_0_t = 0;

K_init_list = [ 1:10, 50, 100, 1000:1000:10000];
K_0_init_list = [ 1:50, 100, 1000, 10000, 100000, 100000];
B = 1;
% % ---------- Check -------------
is_feasible = 0;
for i1 = 1 : length(K_init_list)
    K = K_init_list(i1)*ones(N,1);	T_1 = max(C_n./F_n.*K);	T_2 = max(K);
    for i2 = 1 : length(K_0_init_list)
        K_0 = K_0_init_list(i2);
        
        cons1 = isempty(   find( (C_n./F_n.*K <= T_1)==0, 1 )   );      % Eq. 1
        cons2 = isempty(   find( (K <= T_2*ones(N, 1))==0, 1 )   );      % Eq. 2
        cons3 = isempty(   find( ((B*T_1+C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 <= T_max)==0, 1 )   );
        cons4 = isempty(   find( (c_1/(gamma*K_0*sum(K)) + c_2*gamma^2*T_2^2 ...
            + c_3*gamma/B + c_4*gamma*sum(q_s0sn.*K.^2)/sum(K) <= C_max)==0, 1)   );
        
        if cons1 && cons2 && cons3 && cons4
            is_feasible = 1;    fprintf('Feasible initial point found!\n');   break;
        end
        if is_feasible
            break;
        end
    end
    if is_feasible
        break;
    end
end

if is_feasible == 0
    fprintf('Error: Feasible point not found!  \n');
    pause(1);
    return;
end

% ------------- Algorithm --------------
energy_rcd = [];
norm_rcd = [];
while norm([K_t; K_0_t]-[K; K_0], 2)>0.1
    K_t = K;
    K_0_t = K_0;
    
    beta = K_t./sum(K_t);
    
    cvx_begin gp quiet
    cvx_solver SeDuMi
    variables K(N,1) K_0 T_1 T_2
    minimize (   K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0)  )
    subject to
    C_n./F_n.*K/T_1  <= ones(N,1);
    K/T_2 <= ones(N,1);
    (C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 + B*K_0*T_1 <= T_max;
    c_1/(gamma*K_0*prod((K./beta).^beta)) + c_2*gamma^2*T_2^2 + c_3*gamma/B ...
        + c_4*gamma*sum(q_s0sn.*K.^2)/prod((K./beta).^beta) <= C_max;
%     K>= 0.5*ones(N, 1);
%     K_0 >= 0.5;
%     B >= 0.5;
    cvx_end
    
    energy_rcd = [energy_rcd, K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0)];
    norm_rcd = [norm_rcd, norm([K_t; K_0_t]-[K; K_0], 2)];
end
energy = K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);

% % K(K<1) = 1;
% K_0(K_0<1) = 1;
% 
% K_0_tmp = [ ceil(K_0), floor(K_0) ];   rcv_error = 10000000000;
% K1_tmp = [ ceil( K(1:floor(N/2)) ), floor( K(1:floor(N/2)) ) ]; B_tmp = [ ceil(B), floor(B) ]; K_0_tmp = [ ceil(K_0), floor(K_0) ];   rcv_error = 10000000000;
% K2_tmp = [ ceil( K(floor(N/2)+1:end) ), floor( K(floor(N/2)+1:end) ) ];
% 
% for i1 = 1 : size(K1_tmp, 2)
%     K(1:floor(N/2)) = K1_tmp(:, i1);
%     for i2 = 1 : size(K2_tmp, 2)
%         K(floor(N/2)+1:end) = K2_tmp(:, i2);
%         for j = 1 : length(K_0_tmp)
%             K_0 = K_0_tmp(j);
%             energy_tmp = K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);
%             if abs(energy_tmp-energy_rcd(end)) < rcv_error
%                 energy = energy_tmp;
%                 rcv_error = abs(energy_tmp-energy_rcd(end));
%             end
%         end
%     end
% end
rcv_error = 10000000000;
K_0_tmp = K_0;  K_tmp = K;
for i = 1 : 2 : 1000
    if floor( K_0_tmp ) + i - 500 > 0
        K_0 = floor(K_0_tmp) + i - 500;
    end
    for j1 = 1 : 20
        if floor( K_tmp(1) ) + j1 - 10 >= 0
            K(1:floor(N/2)) = floor( K_tmp(1:floor(N/2)) ) + j1 - 10;
        end
        for j2 = 1 : 20
            if floor( K_tmp(end) ) + j2 - 10 >= 0
                K(floor(N/2)+1:end) = floor( K_tmp(floor(N/2)+1:end) ) + j2 - 10;
            end
            energy_tmp = K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);
            T_1 = max(C_n./F_n.*K);	T_2 = max(K);
            cons3 = isempty(   find( ((B*T_1+C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 <= T_max)==0, 1 )   );
            cons4 = isempty(   find( (c_1/(gamma*K_0*sum(K)) + c_2*gamma^2*T_2^2 ...
                + c_3*gamma/B + c_4*gamma*sum(q_s0sn.*K.^2)/sum(K) <= C_max)==0, 1)   );
            if cons3 && cons4 && energy_tmp-energy_rcd(end) < rcv_error
                rcv_error = energy_tmp-energy_rcd(end);
                energy = energy_tmp;
                fprintf('OK!');
            end
        end
        
    end
end
fprintf('\n');