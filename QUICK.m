function T = QUICK

rho = 997;
k = 0.595;

N = 5;

L = 0.1;

T0 = 150;
Tl = 50;

dx = L/N;

u = 1*ones(1,N);

ue = u;
uw = u;

Fe = rho*ue;
Fw = rho*uw;

De = k/dx;
Dw = k/dx;

a_ee = zeros(1,N);
a_ww = -Fw/8;
a_ww(1) = 0;
a_ww(2) = 0;

a_w = Dw + 6*Fw/8 + 1*Fe/8;
a_w(1) = 0;
a_w(2) = Dw + 7*Fw(2)/8 + 1*Fe(2)/8;
a_w(N) = Dw + 1/3*De + 6/8*Fw(N);

a_e = De - 3*Fe/8;
a_e(N) = 0;
a_e(1) = De + 1*Dw/3 - 3*Fe(1)/8;

S_u = zeros(1,N);
S_u(1) = (8/3*Dw + 2/8*Fe(1) + Fw(1))*T0;
S_u(2) = -0.25*Fw(2)*T0;
S_u(N) = (8/3*De - Fe(N))*Tl;

S_p = zeros(1,N);
S_p(1) = -(8/3*Dw + 2/8*Fe(1) + Fw(1));
S_p(2) = 0.25*Fw(2);
S_p(N) = -(8/3*De - Fe(N));

a_p = a_e + a_w + a_ww + a_ee + (Fe - Fw) - S_p;

A = zeros(N,N);

for i = 1 : N
    if i == 1
        A(i,i) = a_p(i);
        A(i,i+1) = -a_e(i);
        A(i,i+2) = -a_ee(i);
    elseif i == 2
        A(i,i-1) = -a_w(i);
        A(i,i) = a_p(i);
        A(i,i+1) = -a_e(i);
        A(i,i+2) = -a_ee(i);
    elseif i == N
        A(i,i-2) = -a_ww(i);
        A(i,i-1) = -a_w(i);
        A(i,i) = a_p(i);
    elseif i == N-1
        A(i,i-2) = -a_ww(i);
        A(i,i-1) = -a_w(i);
        A(i,i) = a_p(i);
        A(i,i+1) = -a_e(i);
    else
        A(i,i-2) = -a_ww(i);
        A(i,i-1) = -a_w(i);
        A(i,i) = a_p(i);
        A(i,i+1) = -a_e(i);
        A(i,i+2) = -a_ee(i);
    end
end

T = A\(S_u');
T = T';
