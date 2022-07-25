function JacobiVector

% Gauss-Seidel Algorithm for Diffusion 1D Steady-State without
% internal energy generation and given edge temperatures
% Author: Conrado Chiarello
% NUEM - Multiphase Flow Research Center
% Professor: Dr. Paulo H. D. Santos
tic
% Given temperatures
Tp1 = 150;
Tp2 = 50;

% Number of elements
N = 1000;

a_e = ones(1,N);
a_e(end) = 0;
a_w = ones(1,N);
a_w(1) = 0;
S_u = zeros(1,N);
S_u(1) = 2*Tp1;
S_u(end) = 2*Tp2;
S_p = zeros(1,N);
S_p(1) = -2;
S_p(end) = -2;
a_p = a_e + a_w - S_p;

% Initial guess
T = ones(1,N)*(0.5*(Tp1 + Tp2));
T_e = [T(2:N), 0];
T_w = [0, T(1:N-1)];

while (1)
    T = (a_e.*T_e + a_w.*T_w + S_u)./a_p;
    T_e = [T(2:N), 0];
    T_w = [0, T(1:N-1)];
    RMS = sqrt(sum((S_u + a_e.*T_e + a_w.*T_w - a_p.*T).^2));
    if RMS < 1E-3
        break
    end
end
toc
