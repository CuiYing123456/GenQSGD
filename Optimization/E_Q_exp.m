function [ energy] = E_Q_exp( N, s_n, s_0, D, alpha, ...
    C_0, F_0, p_0, r_0, F_n, C_n, p_n, r_n, ...
    Q, gamma, rho_E, G, L, sigma, C_max, T_max )

energy = Inf;

q_0 = min(D/s_0^2, sqrt(D)/s_0);
q_n = min([D./(s_n.^2), sqrt(D)./s_n], [], 2);
M_n = 32 + D + D*log2(s_n);
M_0 = 32 + D + D*log2(s_0);

a_1 = (1-rho_E)/gamma;
a_2 = gamma^2/(1+rho_E+rho_E^2);
a_3 = gamma/(1+rho_E);
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

K_init_list = [ 1:50, 60:10:120, 140:20:800, 900:100:2000];
B_init_list = [ 1:50, 60:10:100, 200:100:1000, 2000:1000:10000];
X_0_init_list = [0.9999, 0.9995, 0.999, 0.995, 0.99, 0.95, 0.9, 0.5, 0.09, 0.05, 0.009, 0.005, 0.0009, 0.0005, 0.00009];

% % ---------- Check -------------
is_feasible = 0;
for i1 = 1 : length(K_init_list)
    K = K_init_list(i1)*ones(N,1);	T_1 = max(C_n./F_n.*K);	T_2 = max(K);
    for i2 = 1 : length(X_0_init_list)
        X_0 = X_0_init_list(i2);
        K_0 = log(X_0)/log(rho_E);
        for i3 = 1 : length(B_init_list)
            B = B_init_list(i3);

            cons1 = isempty(   find( (C_n./F_n.*K <= T_1)==0, 1 )   );      % Eq. 1
            cons2 = isempty(   find( (K <= T_2*ones(N, 1))==0, 1 )   );      % Eq. 2
            cons3 = isempty(   find( ((B*T_1+C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 <= T_max)==0, 1 )   );
            cons4 = isempty(   find( (a_1*c_1/((1-X_0)*sum(K)) + a_2*c_2*(1-X_0^3)*T_2^2 / (1-X_0) ...
                + a_3*(1-X_0^2)/(1-X_0)*(c_3/B+c_4*sum(q_s0sn.*K.^2)/sum(K)) <= C_max)==0, 1)   );
            cons5 = ( abs( log(X_0) - K_0*log(rho_E) ) < 1e-6 );

            if cons1 && cons2 && cons3 && cons4 && cons5
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

if is_feasible == 0
    fprintf('Error: Feasible point not found!  \n');
    pause(1);
    return;
end

% ----------- Algorithm ------------
energy_rcd = [];
norm_rcd = [];

while norm([K_t; K_0_t; B_t]-[K; K_0; B], 2)>0.1
    K_t = K;
    K_0_t = K_0;
    B_t = B;
    T_2_t = T_2;
    X_0_t = X_0;

    lambda_0 = (C_max+a_2*c_2*T_2_t^2*X_0_t^3+a_3*c_3*X_0_t^2/B_t)*sum(K_t) + a_3*c_4*X_0_t^2*sum(q_s0sn.*K_t.^2);
    lambda_1 = C_max*K_t/lambda_0;
    lambda_2 = a_2*c_2*T_2_t^2*X_0_t^3*K_t/lambda_0;
    lambda_3 = a_3*c_3*X_0_t^2*K_t/(B_t*lambda_0);
    lambda_4 = a_3*c_4*q_s0sn*X_0_t^2.*K_t.^2/lambda_0;

    theta_1 = K_0_t*log(1/rho_E) / (K_0_t*log(1/rho_E)+1);
    theta_2 = 1 / (K_0_t*log(1/rho_E)+1);

    cvx_begin gp quiet
    cvx_precision low
    cvx_solver SeDuMi
    variables K(N,1) K_0 B T_1 T_2 X_0
    minimize (   K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0)  )
    subject to
    C_n./F_n.*K/T_1  <= ones(N,1);
    K/T_2 <= ones(N,1);
    X_0 / 0.99999 <= 1;
    ( (C_0/F_0+max(M_n./r_n)+M_0/r_0)*K_0 + B*K_0*T_1 ) / T_max <= 1;

    ( a_1*c_1 + (a_2*c_2*T_2^2+a_3*c_3/B+C_max*X_0)*sum(K) + a_3*c_4*sum(q_s0sn.*K.^2) ) ...
        / prod( ((C_max*K./lambda_1).^lambda_1) .* ((a_2*c_2*T_2^2*X_0^3.*K./lambda_2).^lambda_2) ...
        .* ((a_3*c_3*X_0^2.*K./(B.*lambda_3)).^lambda_3) .* ((a_3*c_4*q_s0sn.*X_0.^2.*K.^2./lambda_4).^lambda_4) ) <= 1;

    ( log(1/X_0_t)*X_0+X_0_t )  /  (  X_0 ...
        * (  K_0*log(1/rho_E) / theta_1 )^theta_1 ...
        * ( 1/theta_2 )^theta_2 ) <= 1;

    ( X_0/X_0_t + K_0*log(1/rho_E) ) / (log(1/X_0_t) + 1) <= 1;

    cvx_end

    energy_rcd = [energy_rcd, K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0)];
    norm_rcd = [norm_rcd, norm([K_t; K_0_t; B_t]-[K; K_0; B], 2)];
end
energy = K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);

% K(K<1) = 1;
B(B<1) = 1;
K_0(K_0<1) = 1;

K1_tmp = [ ceil( K(1:floor(N/2)) ), floor( K(1:floor(N/2)) ) ]; B_tmp = [ ceil(B), floor(B) ]; K_0_tmp = [ ceil(K_0), floor(K_0) ];   rcv_error = 10000000000;
K2_tmp = [ ceil( K(floor(N/2)+1:end) ), floor( K(floor(N/2)+1:end) ) ];

for i1 = 1 : 2
    K(1:floor(N/2)) = K1_tmp(:, i1);
    for i2 = 1 : 2
        K(floor(N/2)+1:end) = K2_tmp(:, i2);
        for j = 1 : length(K_0_tmp)
            K_0 = K_0_tmp(j);
            for k = 1 : length(B_tmp)
                B = B_tmp(k);
                energy_tmp = K_0 * (B*sum(alpha*C_n.*F_n.^2.*K) + alpha*C_0*F_0^2 + sum(p_n.*M_n./r_n) + p_0*M_0/r_0);
                if abs(energy_tmp-energy_rcd(end)) < rcv_error
                    energy = energy_tmp;
                    rcv_error = abs(energy_tmp-energy_rcd(end));
                end
            end
        end
    end
end
