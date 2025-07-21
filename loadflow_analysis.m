clc; clear;

%% Step 1: System Data

% Bus data: [Bus, Type (1=Slack, 2=PV, 3=PQ), V, θ, PG, QG, PL, QL]
bus_data = [
    1, 1, 1.0, 0,    0,  0,   0,  0;
    2, 2, 1.0, 0,  100, 50,  70, 30;
    3, 3, 1.0, 0,    0,  0,  90, 45;
    4, 3, 1.0, 0,    0,  0, 110, 60;
    5, 3, 1.0, 0,    0,  0,  80, 40;
];

% Line data: [From, To, R, X, B]
line_data = [
    1, 2, 0.01, 0.1,  0.02;
    1, 3, 0.02, 0.15, 0.025;
    2, 3, 0.01, 0.1,  0.015;
    2, 4, 0.015,0.12, 0.02;
    3, 5, 0.02, 0.1,  0.03;
];

n_buses = size(bus_data, 1);

%% Step 2: Ybus Formation
Ybus = zeros(n_buses);

for k = 1:size(line_data,1)
    from = line_data(k,1);
    to = line_data(k,2);
    R = line_data(k,3);
    X = line_data(k,4);
    B = line_data(k,5);
    
    Z = R + 1i*X;
    Y = 1/Z;
    
    Ybus(from, from) = Ybus(from, from) + Y + 1i*B;
    Ybus(to, to)     = Ybus(to, to) + Y + 1i*B;
    Ybus(from, to)   = Ybus(from, to) - Y;
    Ybus(to, from)   = Ybus(to, from) - Y;
end

%% Step 3: Initialization
V = bus_data(:,3);                  % Initial voltages
theta = deg2rad(bus_data(:,4));    % Convert angle to radians
P_gen = bus_data(:,5);
Q_gen = bus_data(:,6);
P_load = bus_data(:,7);
Q_load = bus_data(:,8);

P_spec = (P_gen - P_load)/100;     % Per unit
Q_spec = (Q_gen - Q_load)/100;

bus_type = bus_data(:,2);
PQ_idx = find(bus_type == 3);
PV_idx = find(bus_type == 2);
PQ_PV_idx = [PV_idx; PQ_idx];

max_iter = 20;
tol = 1e-6;

%% Step 4: Newton-Raphson Iteration
for iter = 1:max_iter
    P_calc = zeros(n_buses,1);
    Q_calc = zeros(n_buses,1);
    
    for i = 1:n_buses
        for j = 1:n_buses
            G = real(Ybus(i,j));
            B = imag(Ybus(i,j));
            P_calc(i) = P_calc(i) + V(i)*V(j)*(G*cos(theta(i)-theta(j)) + B*sin(theta(i)-theta(j)));
            Q_calc(i) = Q_calc(i) + V(i)*V(j)*(G*sin(theta(i)-theta(j)) - B*cos(theta(i)-theta(j)));
        end
    end
    
    dP = P_spec - P_calc;
    dQ = Q_spec - Q_calc;
    
    mismatch = [dP(PQ_PV_idx); dQ(PQ_idx)];
    
    if max(abs(mismatch)) < tol
        fprintf('Converged in %d iterations.\n\n', iter);
        break;
    end
    
    % Jacobian
    npv = length(PV_idx);
    npq = length(PQ_idx);
    J11 = zeros(npv+npq, npv+npq); % ∂P/∂θ
    J12 = zeros(npv+npq, npq);     % ∂P/∂V
    J21 = zeros(npq, npv+npq);     % ∂Q/∂θ
    J22 = zeros(npq, npq);         % ∂Q/∂V

    for i = 1:(npv+npq)
        m = PQ_PV_idx(i);
        for j = 1:(npv+npq)
            n = PQ_PV_idx(j);
            if m == n
                for k = 1:n_buses
                    if k ~= m
                        G = real(Ybus(m,k));
                        B = imag(Ybus(m,k));
                        J11(i,j) = J11(i,j) + V(m)*V(k)*(-G*sin(theta(m)-theta(k)) + B*cos(theta(m)-theta(k)));
                    end
                end
            else
                G = real(Ybus(m,n));
                B = imag(Ybus(m,n));
                J11(i,j) = V(m)*V(n)*(G*sin(theta(m)-theta(n)) - B*cos(theta(m)-theta(n)));
            end
        end
    end

    for i = 1:(npv+npq)
        m = PQ_PV_idx(i);
        for j = 1:npq
            n = PQ_idx(j);
            if m == n
                for k = 1:n_buses
                    G = real(Ybus(m,k));
                    B = imag(Ybus(m,k));
                    J12(i,j) = J12(i,j) + V(k)*(G*cos(theta(m)-theta(k)) + B*sin(theta(m)-theta(k)));
                end
                J12(i,j) = J12(i,j) + 2*V(m)*real(Ybus(m,m));
            else
                G = real(Ybus(m,n));
                B = imag(Ybus(m,n));
                J12(i,j) = V(m)*(G*cos(theta(m)-theta(n)) + B*sin(theta(m)-theta(n)));
            end
        end
    end

    for i = 1:npq
        m = PQ_idx(i);
        for j = 1:(npv+npq)
            n = PQ_PV_idx(j);
            if m == n
                for k = 1:n_buses
                    if k ~= m
                        G = real(Ybus(m,k));
                        B = imag(Ybus(m,k));
                        J21(i,j) = J21(i,j) + V(m)*V(k)*(G*cos(theta(m)-theta(k)) + B*sin(theta(m)-theta(k)));
                    end
                end
            else
                G = real(Ybus(m,n));
                B = imag(Ybus(m,n));
                J21(i,j) = -V(m)*V(n)*(G*cos(theta(m)-theta(n)) + B*sin(theta(m)-theta(n)));
            end
        end
    end

    for i = 1:npq
        m = PQ_idx(i);
        for j = 1:npq
            n = PQ_idx(j);
            if m == n
                for k = 1:n_buses
                    G = real(Ybus(m,k));
                    B = imag(Ybus(m,k));
                    J22(i,j) = J22(i,j) + V(k)*(G*sin(theta(m)-theta(k)) - B*cos(theta(m)-theta(k)));
                end
                J22(i,j) = J22(i,j) - 2*V(m)*imag(Ybus(m,m));
            else
                G = real(Ybus(m,n));
                B = imag(Ybus(m,n));
                J22(i,j) = V(m)*(G*sin(theta(m)-theta(n)) - B*cos(theta(m)-theta(n)));
            end
        end
    end

    % Combine Jacobian
    J = [J11 J12; J21 J22];

    % Solve
    delta = J \ mismatch;

    % Update
    theta(PQ_PV_idx) = theta(PQ_PV_idx) + delta(1:(npv+npq));
    V(PQ_idx) = V(PQ_idx) + delta((npv+npq+1):end);
end

%% Step 5: Final Results
fprintf('Bus\t| V (pu)\t| Angle (deg)\n');
disp('----------------------------------------')
for i = 1:n_buses
    fprintf('%d\t| %.4f\t| %.2f\n', i, V(i), rad2deg(theta(i)));
end

%% Step 6: Plotting
figure;
subplot(2,1,1);
plot(1:n_buses, V, '-o');
title('Voltage Magnitudes');
xlabel('Bus'); ylabel('V (pu)');

subplot(2,1,2);
plot(1:n_buses, rad2deg(theta), '-o');
title('Voltage Angles');
xlabel('Bus'); ylabel('Angle (deg)');
