function T = Upwind

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

if u < 0
    a_e = De + Fe;
    a_e(N) = 0;
    
    a_w = Dw*ones(1,N);
    a_w(1) = 0;
    
    S_u = zeros(1,N);
    S_u(1) = 2*Dw*T0;
    S_u(N) = (2*De + Fe(N))*Tl;
    
    S_p = zeros(1,N);
    S_p(1) = -2*Dw;
    S_p(N) = -(2*De + Fe(N));
elseif u > 0
    a_e = ones(1,N)*De;
    a_e(N) = 0;
    
    a_w = Dw + Fw;
    a_w(1) = 0;
    
    S_u = zeros(1,N);
    S_u(1) = (2*Dw + Fw(1))*T0;
    S_u(N) = 2*De*Tl;
    
    S_p = zeros(1,N);
    S_p(1) = -(2*Dw + Fw(1));
    S_p(N) = -2*De;
end

a_p = a_e + a_w + (Fe - Fw) - S_p;

T = TDMASolver(a_w, a_e, a_p, S_u);
