function [ energy, K, B, K_0, gamma] = E_Q_opt_stepsize( N, s_n, s_0, D, alpha, ...
    C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
    Q, G, L, sigma, C_max, T_max )

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
B_t = 0;
gamma_t = 0;

K_init_list = [ 1:20, 100:100:3000];
K_0_init_list = [ 1:20, 100:200:5000, 10000, 100000, 100000];
B_init_list = [ 1:20, 100, 1000, 10000, 100000];
gamma_init_list = [0.005:0.005:0.06];
rcv_error = 10000000000;
energy_min = 100000000;
% % ---------- Check -------------
is_feasible = 0;
for i1 = 1 : length(K_init_list)
    K = K_init_list(i1)*ones(N,1);	T_1 = max(C_n./F_n.*K);	T_2 = max(K);
    for i2 = 1 : length(K_0_init_list)
        K_0 = K_0_init_list(i2);
        for i3 = 1 : length(B_init_list)
            B = B_init_list(i3);
            for i4 = 1 : length(gamma_init_list)
                gamma = gamma_init_list(i4);

                cons1 = isempty(   find( (C_n./F_n.*K <= T_1)==0, 1 )   );
                cons2 = isempty(   find( (K <= T_2*ones(N, 1))==0, 1 )   );
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
        if is_feasible
            break;
        end
    end

    if is_feasible == 0
        fprintf('Error: Feasible point not found!  \n');
        pause(1);
        return;
    end
end
% ------------- Algorithm --------------
energy_rcd = [];
norm_rcd = [];
while norm([K_t; K_0_t; B_t; gamma_t]-[K; K_0; B; gamma], 2)>0.01
    K_t = K;
    K_0_t = K_0;
    B_t = B;
    gamma_t = gamma;

    beta = K_t./sum(K_t);

    cvx_begin gp quiet
    cvx_solver SeDuMi
    variables K(N,1) K_0 B gamma T_1 T_2
    minimize (   K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0)  )
    subject to
    C_n./F_n.*K/T_1  <= ones(N,1);
    K/T_2 <= ones(N,1);
    (C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 + B*K_0*T_1 <= T_max;
    c_1/(gamma*K_0*prod((K./beta).^beta)) + c_2*gamma^2*T_2^2 + c_3*gamma/B ...
        + c_4*gamma*sum(q_s0sn.*K.^2)/prod((K./beta).^beta) <= C_max;
    gamma <= 1/L;
    %     K>= 0.5*ones(N, 1);
    %     K_0 >= 0.5;
    %     B >= 0.5;
    cvx_end

    energy_rcd = [energy_rcd, K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0)];
    norm_rcd = [norm_rcd, norm([K_t; K_0_t; B_t; gamma_t]-[K; K_0; B; gamma], 2)];
end
energy_min = min([K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0), energy_min]);

K_0_tmp = K_0;  K_tmp = K;  B_tmp = B;  gamma_tmp = gamma;
for i = 1 : 2 : 2000
    if floor( K_0_tmp ) + i - 1000 > 0
        K_0 = floor(K_0_tmp) + i - 1000;
    end
    for j1 = 1 : 20
        if floor( K_tmp(1) ) + j1 - 10 >= 0
            K(1:floor(N/2)) = floor( K_tmp(1:floor(N/2)) ) + j1 - 10;
        end
        for j2 = 1 : 20
            if floor( K_tmp(end) ) + j2 - 10 >= 0
                K(floor(N/2)+1:end) = floor( K_tmp(floor(N/2)+1:end) ) + j2 - 10;
            end
            for k = 1 : 20
                if floor(B_tmp) + k - 10 > 0
                    B = floor(B_tmp) + k - 10;
                end
                for n = 1 : 20
                    if gamma_tmp - 0.02 + 0.02*n/10 >= 0
                        gamma = gamma_tmp - 0.02 + 0.02*n/10;
                    end
                    energy_tmp = K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);
                    T_1 = max(C_n./F_n.*K);	T_2 = max(K);
                    cons3 = isempty(   find( ((B*T_1+C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 <= T_max)==0, 1 )   );
                    cons4 = isempty(   find( (c_1/(gamma*K_0*sum(K)) + c_2*gamma^2*T_2^2 ...
                        + c_3*gamma/B + c_4*gamma*sum(q_s0sn.*K.^2)/sum(K) <= C_max)==0, 1)   );
                    if cons3 && cons4 && energy_tmp-energy_rcd(end) < rcv_error
                        energy = energy_tmp;
                        rcv_error = energy_tmp-energy_rcd(end);
                        fprintf('OK!');
                    end
                    if B == 5 && K(1) == 1 && K(end) == 1 && K_0 == 995
                        a=1;
                    end
                end
            end
        end
    end
end
fprintf('\n');
end