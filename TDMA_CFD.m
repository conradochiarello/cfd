function TDMA_CFD(N)

% TriDiagonal Matrix Algorithm for Diffusion 1D Steady-State without
% internal energy generation and given edge temperatures
% Author: Conrado Chiarello
% NUEM - Multiphase Flow Research Center
% Professor: Dr. Paulo H. D. Santos

tic
% Given temperatures
Tp1 = 150;
Tp2 = 50;

% Initializing variables
a_p = zeros(N,1);
a_e = zeros(N,1);
a_w = zeros(N,1);
S_u = zeros(N,1);
S_p = zeros(N,1);
Q = zeros(N,1);
P = zeros(N,1);
A = zeros(N,1);
B = zeros(N,1);
C = zeros(N,1);
D = zeros(N,1);

% Main loop
for i = 1 : N
   if i == 1
       a_e(i) = 1;
       a_w(i) = 0;
       S_u(i) = 2*Tp1;
       S_p(i) = -2;
       a_p(i) = a_e(i) + a_w(i) - S_p(i);

       A(i) = a_p(i);
       B(i) = -a_e(i);
       C(i) = -a_w(i);
       D(i) = S_u(i);
       
       P(i) = -B(i)/A(i);
       Q(i) = D(i)/A(i);
   elseif i < N
       a_e(i) = 1;
       a_w(i) = 1;
       S_u(i) = 0;
       S_p(i) = 0;
       a_p(i) = a_e(i) + a_w(i) - S_p(i);
       
       A(i) = a_p(i);
       B(i) = -a_e(i);
       C(i) = -a_w(i);
       D(i) = S_u(i);
       
       P(i) = -B(i)/(A(i) + C(i)*P(i-1));
       Q(i) = (D(i) - C(i)*Q(i-1))/(A(i) + C(i)*P(i-1));
   elseif i == N
       a_e(i) = 0;
       a_w(i) = 1;
       S_u(i) = 2*Tp2;
       S_p(i) = -2;
       a_p(i) = a_e(i) + a_w(i) - S_p(i);
       
       A(i) = a_p(i);
       B(i) = -a_e(i);
       C(i) = -a_w(i);
       D(i) = S_u(i);
       
       P(i) = -B(i)/(A(i) + C(i)*P(i-1));
       Q(i) = (D(i) - C(i)*Q(i-1))/(A(i) + C(i)*P(i-1));
   end  
end

T(N) = Q(N);

for i = N : -1 : 2
        T(i-1) = P(i-1)*T(i) + Q(i-1);
end
toc
