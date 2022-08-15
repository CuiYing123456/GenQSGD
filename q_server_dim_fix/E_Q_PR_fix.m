function [ energy] = E_Q_PR_fix( N, s_n, s_0, D, alpha, ...
    C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
    Q, gamma, rho_D, G, L, sigma, C_max, T_max )

energy = Inf;

q_0 = min(D/s_0^2, sqrt(D)/s_0);
q_n = min([D./(s_n.^2), sqrt(D)./s_n], [], 2);
M_n = 32 + D + D*log2(s_n);
M_0 = 32 + D + D*log2(s_0);

b_1 = 1/(rho_D*gamma);
b_2 = rho_D^2*gamma^2/(rho_D+1)^3 + rho_D^2*gamma^2/(2*(rho_D+1)^2);
b_3 = rho_D*gamma/(rho_D+1)^2 + rho_D*gamma/(rho_D+1);
c_1 = 2*N*Q;
c_2 = 4*G^2*L^2;
c_3 = L*sigma^2/N;
c_4 = 2*L*G^2;
q_s0sn = q_0+q_n+q_0*q_n;

%% =============== Solution: Proposed =================
K_init_list = [ 1:10, 20:10:100, 200:500:2000];
K_0_init_list = [ 1:10, 20:10:100, 300:500:3000];
B = 1;
for i1 = 1 : length(K_init_list)
    K = K_init_list(i1)*ones(N,1);	T_1 = max(C_n./F_n.*K);	T_2 = max(K);
    for i2 = 1 : length(K_0_init_list)
        K_0 = K_0_init_list(i2);
        cons3 = isempty(   find( ((B*T_1+C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 <= T_max)==0, 1 )   );
        cons4 = ( b_1*c_1/(sum(K)) + b_2*c_2*T_2^2 + b_3*c_3/B + b_3*c_4*sum(q_s0sn.*K.^2)/sum(K) <= C_max*log( (K_0+rho_D+1) / (rho_D+1) ) );
        
        if cons3 && cons4
            fprintf('Feasible!\n');
            energy_tmp = K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);
            if energy_tmp < energy
                energy = energy_tmp;
            end
        end
    end
end
