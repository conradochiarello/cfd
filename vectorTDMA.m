function vectorTDMA
tic
N = 500000000;

Tp1 = 150;
Tp2 = 50;

P = zeros(N,1);
Q = zeros(N,1);

a_w = ones(1,N);
a_w(1) = 0;

a_e = ones(1,N);
a_e(N) = 0;

S_u = zeros(1,N);
S_u(1) = 2*Tp1;
S_u(N) = 2*Tp2;

S_p = zeros(1,N);
S_p(1) = -2;
S_p(N) = -2;

a_p = a_e + a_w - S_p;

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
toc
