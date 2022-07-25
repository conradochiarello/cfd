function [u, p] = SIMPLE

% Fluid properties
rho = 1;

% Number of points
N = 5;

% Tolerance
tol = 1E-4;

% Relaxation coefficient
wv = 0.5;
wp = 1/N;

% Areas (pressure and velocity)
A = linspace(0.5, 0.1, N);
A_v = linspace(0.45, 0.15, N-1);

% Estimated mass flow
mdot = 1;

% Given pressures
p0 = 10;
pend = 0;

% Estimated velocity and pressure fields
u = mdot./(rho*A_v);
p = linspace(p0, pend, N);

% Preallocating variables for code performance
ue = zeros(1,N-1);
uw = zeros(1,N-1);
Fe = zeros(1,N-1);
Fw = zeros(1,N-1);
a_e = zeros(1,N-1);
a_w = zeros(1,N-1);
a_p = zeros(1,N-1);
S_u = zeros(1,N-1);
ustar = zeros(1,N-1);
bprime = zeros(1,N-1);
awp = zeros(1,N-2);
aep = zeros(1,N-2);
app = zeros(1,N-2);
Fwp = zeros(1,N-2);
Fep = zeros(1,N-2);
ucalc = zeros(1,N-1);
pCorrection = zeros(1,N);

while(true)
    for i = 1 : N-1
        if i == 1
            Fw(i) = rho*u(i)*A_v(i);
            Fe(i) = rho*0.5*(u(i+1) + u(i))*A(i+1);
            a_p(i) = Fe(i) + Fw(i)*0.5*(A_v(i)/A(i))^2;
            S_u(i) = (p(i) - p(i+1))*A_v(i) + Fw(i)*A_v(i)/A(i)*u(i);
        elseif i == N - 1
            Fw(i) = rho*0.5*(u(i) + u(i-1))*A(i);
            Fe(i) = rho*u(i)*A_v(i);
            a_w(i) = Fw(i);
            a_p(i) = a_w(i) + a_e(i) + Fe(i) - Fw(i);
            S_u(i) = (p(i) - p(i+1))*A_v(i);
        else
            ue(i) = 0.5*(u(i) + u(i+1));
            uw(i) = 0.5*(u(i-1) + u(i));
            Fe(i) = rho*ue(i)*A(i+1);
            Fw(i) = rho*uw(i)*A(i);
            a_w(i) = Fw(i);
            a_p(i) = a_w(i) + (Fe(i) - Fw(i));
            S_u(i) = (p(i) - p(i+1))*A_v(i);
        end
    end
    
    ustar(1) = S_u(1)/a_p(1);
    
    for i = 2 : N-1
        ustar(i) = (a_w(i)*ustar(i-1) + S_u(i))/a_p(i);
    end
    
    d = A_v./a_p;
    
    
    for i = 1 : N-2
        awp(i) = rho*d(i)*A_v(i);
        aep(i) = rho*d(i+1)*A_v(i+1);
        Fwp(i) = rho*ustar(i)*A_v(i);
        Fep(i) = rho*ustar(i+1)*A_v(i+1);
        
        app(i) = awp(i) + aep(i);
        bprime(i) = Fwp(i) - Fep(i);
    end
    
    pCor = TDMASolver(awp, aep, app, bprime);
    for i = 1 : N
        if i == 1
            pCorrection(i) = 0;
        elseif i == N
            pCorrection(i) = 0;
        else
            pCorrection(i) = pCor(i-1);
        end
    end
    
    p = p + wp*pCorrection;
    
    for i = 1 : N-1
        ucalc(i) = ustar(i) + d(i)*(pCorrection(i) - pCorrection(i+1));
    end
    
    p(1) = p0 - 0.5*rho*ucalc(1)^2*(A_v(1)/A(1))^2;
    u = (1-wv)*u + wv*ucalc;
    
    for i = 1 : N-1
        if i == 1
            Fw(i) = rho*u(i)*A_v(i);
            Fe(i) = rho*0.5*(u(i+1) + u(i))*A(i+1);
            a_p(i) = Fe(i) + Fw(i)*0.5*(A_v(i)/A(i))^2;
            S_u(i) = (p(i) - p(i+1))*A_v(i) + Fw(i)*A_v(i)/A(i)*u(i);
        elseif i == N - 1
            Fw(i) = rho*0.5*(u(i) + u(i-1))*A(i);
            Fe(i) = rho*u(i)*A_v(i);
            a_w(i) = Fw(i);
            a_p(i) = a_w(i) + a_e(i) + Fe(i) - Fw(i);
            S_u(i) = (p(i) - p(i+1))*A_v(i);
        else
            ue(i) = 0.5*(u(i) + u(i+1));
            uw(i) = 0.5*(u(i-1) + u(i));
            Fe(i) = rho*ue(i)*A(i+1);
            Fw(i) = rho*uw(i)*A(i);
            a_w(i) = Fw(i);
            a_p(i) = a_w(i) + (Fe(i) - Fw(i));
            S_u(i) = (p(i) - p(i+1))*A_v(i);
        end
    end
    
    ustar(1) = S_u(1)/a_p(1);
    
    for i = 2 : N-1
        ustar(i) = (a_w(i)*ustar(i-1) + S_u(i))/a_p(i);
    end
    
    d = A_v./a_p;
    
    
    for i = 1 : N-2
        awp(i) = rho*d(i)*A_v(i);
        aep(i) = rho*d(i+1)*A_v(i+1);
        Fwp(i) = rho*ustar(i)*A_v(i);
        Fep(i) = rho*ustar(i+1)*A_v(i+1);
        
        app(i) = awp(i) + aep(i);
        bprime(i) = Fwp(i) - Fep(i);
    end
    
    if max(abs(bprime)) < tol
        break
    end
end
