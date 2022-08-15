function [ energy ] = E_Q_PM( N, s_n, s_0, D, alpha, ...
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

%% =============== Solution: PM-SGD =================
% ------------- Initialize --------------
K_0_t = 0;
B_t = 0;

K_0_init_list = [ 1, 2, 10, 400, 4000, 40000, 100000];
B_init_list = [ 1, 2, 10, 100, 1000, 10000];
% % ---------- Check -------------
is_feasible = 0;
K = 1;

for i2 = 1 : length(K_0_init_list)
    K_0 = K_0_init_list(i2);
    for i3 = 1 : length(B_init_list)
        B = B_init_list(i3);
        cons3 = isempty(   find( ((B*K*max(C_n./F_n)+C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 <= T_max)==0, 1 )   );
        cons4 = ( b_1*c_1/(N*K) + b_2*c_2*K^2 + b_3*c_3/B + b_3*c_4*sum(q_s0sn)*K/N <= C_max*log( (K_0+rho_D+1) / (rho_D+1) ) );
        
        if cons3 && cons4
            is_feasible = 1;    fprintf('Feasible initial point found!\n');   break;
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
while norm([K_0_t; B_t]-[K_0; B], 2)>0.01
    K_0_t = K_0;
    B_t = B;
    
    cvx_begin gp quiet
    cvx_solver SeDuMi
    variables K_0 B
    minimize (   K_0 * (B*K*sum(alpha*C_n.*F_n.^2) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0)  )
    subject to
    (C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 + B*K_0*K*max(C_n./F_n) <= T_max;
    
    (  b_1*c_1/(N*K)  +  b_2*c_2*K^2  +  b_3*c_3/B  +  b_3*c_4*sum(q_s0sn)*K/N ...
         +  C_max*K_0_t^2/(K_0*(K_0_t+rho_D+1))  ) / (  C_max*( log((K_0_t+rho_D+1)/(rho_D+1)) + ( K_0_t/(K_0_t+rho_D+1) ) )  ) <= 1;
     
    cvx_end
    
    energy_rcd = [energy_rcd, K_0 * (B*K*sum(alpha*C_n.*F_n.^2) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0) ];
    norm_rcd = [norm_rcd, norm([K_0_t; B_t]-[K_0; B], 2)];
end
energy = K_0 * (B*K*sum(alpha*C_n.*F_n.^2) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);

% B(B<1) = 1;
% K_0(K_0<1) = 1;
% 
% B_tmp = [ ceil(B), floor(B) ];  K_0_tmp = [ ceil(K_0), floor(K_0) ];   rcv_error = 10000000000;
% 
% for i = 1 : length(B_tmp)
%     B = B_tmp(i);
%     for j = 1 : length(K_0_tmp)
%         K_0 = K_0_tmp(j);
%         energy_tmp = K_0 * (B*K*sum(alpha*C_n.*F_n.^2) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);
%         if abs(energy_tmp-energy_rcd(end)) < rcv_error
%             energy = energy_tmp;
%             rcv_error = abs(energy_tmp-energy_rcd(end));
%         end
%     end
% end
rcv_error = 10000000000;
K_0_tmp = K_0;  B_tmp = B;
for i = 1 : 2 : 1000
    if floor( K_0_tmp ) + i - 500 > 0
        K_0 = floor(K_0_tmp) + i - 500;
    end
    
    for k = 1 : 9
        if floor(B_tmp) + k - 5 > 0
            B = floor(B_tmp) + k - 5;
        end
        energy_tmp = K_0 * (B*K*sum(alpha*C_n.*F_n.^2) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);
        cons3 = isempty(   find( ((B*K*max(C_n./F_n)+C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 <= T_max)==0, 1 )   );
        cons4 = ( b_1*c_1/(N*K) + b_2*c_2*K^2 + b_3*c_3/B + b_3*c_4*sum(q_s0sn)*K/N <= C_max*log( (K_0+rho_D+1) / (rho_D+1) ) );
        if cons3 && cons4 && energy_tmp-energy_rcd(end) < rcv_error
            rcv_error = energy_tmp-energy_rcd(end);
            energy = energy_tmp;
            fprintf('OK!');
        end
    end
    
end
fprintf('\n');
