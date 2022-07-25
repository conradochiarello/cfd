function T = TDMASolver(a_w, a_e, a_p, S_u)

% TriDiagonal Matrix Algorithm
% Author: Conrado Chiarello
% NUEM - Multiphase Flow Research Center
% Professor: Dr. Paulo H. D. Santos

N = length(a_p);

P = zeros(N,1);
Q = zeros(N,1);

A = a_p;
B = -a_e;
C = -a_w;
D = S_u;

P(1) = -B(1)/A(1);
Q(1) = D(1)/A(1);

for i = 2 : N
    P(i) = -B(i)/(A(i) + C(i)*P(i-1));
    Q(i) = (D(i) - C(i)*Q(i-1))/(A(i) + C(i)*P(i-1));
end

T(N) = Q(N);

for i = N : -1 : 2
        T(i-1) = P(i-1)*T(i) + Q(i-1);
end
