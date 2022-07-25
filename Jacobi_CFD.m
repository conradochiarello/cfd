function T = Jacobi_CFD(N)

% Jacobi Algorithm for Diffusion 1D Steady-State without
% internal energy generation and given edge temperatures
% Author: Conrado Chiarello
% NUEM - Multiphase Flow Research Center
% Professor: Dr. Paulo H. D. Santos

tic

% Given temperatures
Tp1 = 150;
Tp2 = 50;

% Initializing variables for code performance
a_p = zeros(N,1);
a_e = zeros(N,1);
a_w = zeros(N,1);
S_u = zeros(N,1);
S_p = zeros(N,1);
res = zeros(N,1);
T = zeros(N,1);

% Initial guess
T_old = ones(N,1)*(0.5*(Tp1 + Tp2));

while (1)
    for i = 1 : N
        if i == 1
            a_e(i) = 1;
            a_w(i) = 0;
            S_u(i) = 2*Tp1;
            S_p(i) = -2;
            a_p(i) = a_e(i) + a_w(i) - S_p(i);
            
            T(i) = (a_e(i)*T_old(i+1) + S_u(i))/a_p(i);
            
        elseif i < N
            a_e(i) = 1;
            a_w(i) = 1;
            S_u(i) = 0;
            S_p(i) = 0;
            a_p(i) = a_e(i) + a_w(i) - S_p(i);
            
            T(i) = (a_e(i).*T_old(i+1) + a_w(i).*T_old(i-1) + S_u(i))./a_p(i);
            
        elseif i == N
            a_e(i) = 0;
            a_w(i) = 1;
            S_u(i) = 2*Tp2;
            S_p(i) = -2;
            a_p(i) = a_e(i) + a_w(i) - S_p(i);
            
            T(i) = (a_w(i)*T_old(i-1) + S_u(i))/a_p(i);
        end
    end
    for i = 1 : N
        if i == 1
            res(i) = S_u(i) + a_e(i)*T(i+1) - a_p(i)*T(i);
        elseif i < N
            res(i) = S_u(i) + a_e(i)*T(i+1) + a_w(i)*T(i-1) - a_p(i)*T(i);
        elseif i == N
            res(i) = S_u(i) + a_w(i)*T(i-1) - a_p(i)*T(i);
        end   
    end
    RMS = sqrt(sum(res.^2));
    if RMS < 1E-3
        break
    end
    T_old = T;
end

toc
