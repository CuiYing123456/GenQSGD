function [ energy ] = E_Q_PM_fix( N, s_n, s_0, D, alpha, ...
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

%% =============== Solution: PM-SGD-fix =================
K = 1;
K_0_init_list = [ 1:20, 20:50:200, 200:600:3000];
B_init_list = [ 1:10, 10:50:200, 200:100:1000];
for i2 = 1 : length(K_0_init_list)
    K_0 = K_0_init_list(i2);
    for i3 = 1 : length(B_init_list)
        B = B_init_list(i3);
        cons3 = isempty(   find( ((B*K*max(C_n./F_n)+C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 <= T_max)==0, 1 )   );
        cons4 = isempty(   find( (c_1/(gamma*K_0*N*K) + c_2*gamma^2*K^2 ...
            + c_3*gamma/B + c_4*gamma*sum(q_s0sn)*K/N <= C_max)==0, 1)   );
        
        if cons3 && cons4
            fprintf('Feasible!\n');
            energy_tmp = K_0 * (B*K*sum(alpha*C_n.*F_n.^2) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);
            if energy_tmp < energy
                energy = energy_tmp;
            end
        end
    end
end