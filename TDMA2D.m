function T = TDMA2D(N)

% Given initial wall temperatures
Tp1 = 25; % Tw
Tp2 = 25; % Te
Tp3 = 25; % Ts
Tp4 = 50; % Tn

% Since k = constant and L = H, they're omitted for this code

tolerance = 1E-3;

% Number of nodes in each row/column is an input

% Initializing P and Q for code performance
P = zeros(N,N);
Q = zeros(N,N);

% Initial guess
T_ini = 0.25*(Tp1 + Tp2 + Tp3 + Tp4);
T = ones(N)*T_ini;
T_n = [zeros(1,N); T_ini*ones(N-1, N)];
T_s = [T_ini*ones(N-1, N); zeros(1,N)];

a_s = [zeros(1, N); ones(N-1, N)];
a_n = [ones(N-1, N); zeros(1, N)];
a_w = [zeros(N, 1), ones(N, N-1)];
a_e = [ones(N, N-1), zeros(N, 1)];
S_u = 2*([ones(N,1)*Tp1, zeros(N, N-1)] + [zeros(N-1, N); Tp3*ones(1, N)] + [Tp4*ones(1,N); zeros(N-1, N)] + [zeros(N, N-1), Tp2*ones(N,1)]);
S_p = -2*([ones(N,1), zeros(N, N-1)] + [zeros(N-1, N); ones(1, N)] + [ones(1,N); zeros(N-1, N)] + [zeros(N, N-1), ones(N,1)]);
a_p = a_e + a_w + a_n + a_s - S_p;

% Main loop, stops at error less than tolerance
while(1)
    
    % First sweep: base to the top
    A = a_p;
    B = -a_e;
    C = -a_w;
    D = a_n.*T_n + a_s.*T_s + S_u;
    P(:,1) = -B(:,1)./A(:,1);
    Q(:,1) = D(:,1)./A(:,1);
    
    % Minor loop to evaluate P(i,j) and Q(i,j)
    for i = 1 : N
        for j = 2 : N
            P(i,j) = -B(i,j)/(A(i,j) + C(i,j)*P(i,j-1));
            Q(i,j) = (D(i,j) - C(i,j)*Q(i, j-1))/(A(i,j) + C(i,j)*P(i,j-1));
        end
    end
    
    T(:,N) = Q(:,N);
    
    % Minor loop to evaluate T(i,j)
    for i = N : -1 : 1
        for j = N : -1 : 2
            T(i, j-1) = P(i, j-1)*T(i,j) + Q(i,j-1);
        end
    end
    
    % Update T_e, T_w, T_n and T_s
    T_s = [zeros(1,N); T(1:N-1,1:N)];
    T_n = [T(2:N,1:N); zeros(1,N)];
    
    % End of the first sweep
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Second sweep: top to the base
    D = a_n.*T_n + a_s.*T_s + S_u;
    P(:,1) = -B(:,1)./A(:,1);
    Q(:,1) = D(:,1)./A(:,1);
    
    % Minor loop to evaluate P(i,j) and Q(i,j)
    for i = 1 : N
        for j = 2 : N
            P(i,j) = -B(i,j)/(A(i,j) + C(i,j)*P(i,j-1));
            Q(i,j) = (D(i,j) - C(i,j)*Q(i, j-1))/(A(i,j) + C(i,j)*P(i,j-1));
        end
    end
    
    T(:,N) = Q(:,N);
    
    % Minor loop to evaluate T(i,j)
    for i = N : -1 : 1
        for j = N : -1 : 2
            T(i, j-1) = P(i, j-1)*T(i,j) + Q(i,j-1);
        end
    end
    
    % Update T_e, T_w, T_n and T_s
    T_w = [zeros(N,1), T(1:N, 1:N-1)];
    T_e = [T(1:N, 2:N), zeros(N,1)];
    
    % End of the second sweep
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Third sweep: left to the right
    B = -a_n;
    C = -a_s;
    D = a_w.*T_w + a_e.*T_e + S_u;
    P(1,:) = -B(1,:)./A(1,:);
    Q(1,:) = D(1,:)./A(1,:);
    
    % Minor loop to evaluate P(i,j) and Q(i,j)
    for j = 1 : N
        for i = 2 : N
            P(i,j) = -B(i,j)/(A(i,j) + C(i,j)*P(i-1,j));
            Q(i,j) = (D(i,j) - C(i,j)*Q(i-1, j))/(A(i,j) + C(i,j)*P(i-1,j));
        end
    end
    
    T(N,:) = Q(N,:);
    
    % Minor loop to evaluate T(i,j)
    for j = N : -1 : 1
        for i = N : -1 : 2
            T(i-1, j) = P(i-1, j)*T(i,j) + Q(i-1,j);
        end
    end
    
    % Update T_e, T_w, T_n and T_s
    T_w = [zeros(N,1), T(1:N, 1:N-1)];
    T_e = [T(1:N, 2:N), zeros(N,1)];
    
    % End of the third sweep
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Fourth sweep: right to the left
    D = a_w.*T_w + a_e.*T_e + S_u;
    P(1,:) = -B(1,:)./A(1,:);
    Q(1,:) = D(1,:)./A(1,:);
    
    % Minor loop to evaluate P(i,j) and Q(i,j)
    for j = 1 : N
        for i = 2 : N
            P(i,j) = -B(i,j)/(A(i,j) + C(i,j)*P(i-1,j));
            Q(i,j) = (D(i,j) - C(i,j)*Q(i-1, j))./(A(i,j) + C(i,j)*P(i-1,j));
        end
    end
    
    T(N,:) = Q(N,:);
    
    % Minor loop to evaluate T(i,j)
    for j = N : -1 : 1
        for i = N : -1 : 2
            T(i-1, j) = P(i-1, j)*T(i,j) + Q(i-1,j);
        end
    end
    
    % Update T_e, T_w, T_n and T_s
    T_w = [zeros(N,1), T(1:N, 1:N-1)];
    T_e = [T(1:N, 2:N), zeros(N,1)];
    T_s = [zeros(1,N); T(1:N-1,1:N)];
    T_n = [T(2:N,1:N); zeros(1,N)];
    
    % Evaluate error
    err = sqrt(sum(sum((S_u + a_e.*T_e + a_w.*T_w + a_n.*T_n + a_s.*T_s - a_p.*T).^2)));
    
    if err < tolerance
        break
    end
    
    % End of the fourth sweep
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

% % Analytic solution for y = H/2 = 0.5
% L = 1;
% m = 100;
% T_analytic = zeros(m,m);
% x = linspace(0,1,m);
% 
% for n = 1 : m
%    C_n = -50*(-1 + cos(n*pi))/(sinh(n*pi)*n*pi);
%    T_analytic = T_analytic + C_n*sin(n*pi*x/L).*sinh(n*pi*0.5/L);
% end
% T_analytic = T_analytic + Tp1;
% figure (1)
% surf(linspace(0,1,N), linspace(0,1,N), T, 'EdgeColor', 'none','DisplayName', 'T(x,y)','FaceColor', 'interp')
% figure (2)
% hold on
% plot(linspace(0,1,N), T(N/2,:))
% plot(x, T_analytic(m/2,:));
% hold off
