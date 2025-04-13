% Load Flow Analysis using Newton-Raphson Method
clc; clear;

% Step 1: Define system parameters
% Bus data: [Bus No., Type (1=Slack, 2=PV, 3=PQ), V (initial), Angle (initial), PG, QG, PL, QL]
bus_data = [
    1, 1, 1.0, 0, 0, 0, 0, 0;  % Slack bus
    2, 2, 1.0, 0, 100, 50, 70, 30; % PV bus
    3, 3, 1.0, 0, 0, 0, 90, 45;  % PQ bus
    4, 3, 1.0, 0, 0, 0, 110, 60; % PQ bus
    5, 3, 1.0, 0, 0, 0, 80, 40;  % PQ bus
];

% Line data: [From Bus, To Bus, R, X, B]
line_data = [
    1, 2, 0.01, 0.1, 0.02;
    1, 3, 0.02, 0.15, 0.025;
    2, 3, 0.01, 0.1, 0.015;
    2, 4, 0.015, 0.12, 0.02;
    3, 5, 0.02, 0.1, 0.03;
];

% Step 2: Initialize Ybus (Admittance Matrix)
n_buses = size(bus_data, 1);
Ybus = zeros(n_buses);

% Populate Ybus matrix based on line data
for i = 1:size(line_data, 1)
    from = line_data(i, 1);
    to = line_data(i, 2);
    R = line_data(i, 3);
    X = line_data(i, 4);
    B = line_data(i, 5);
    
    % Admittance of the line
    Z = R + 1i * X;  % Impedance
    Y = 1 / Z;  % Admittance
    
    % Off-diagonal elements of Ybus
    Ybus(from, to) = Ybus(from, to) - Y;
    Ybus(to, from) = Ybus(to, from) - Y;
    
    % Diagonal elements of Ybus
    Ybus(from, from) = Ybus(from, from) + Y + 1i * B;
    Ybus(to, to) = Ybus(to, to) + Y + 1i * B;
end

% Step 3: Initialize Voltage Magnitudes and Angles
V = bus_data(:, 3);   % Initial voltage magnitudes
theta = bus_data(:, 4);  % Initial voltage angles in radians
V_pu = V;  % Store per unit voltage

% Step 4: Iterative solution using Newton-Raphson method
max_iter = 20;  % Maximum number of iterations
tolerance = 1e-6;  % Convergence tolerance
converged = false;

for iter = 1:max_iter
    % Calculate power mismatch
    P_mismatch = zeros(n_buses, 1);
    Q_mismatch = zeros(n_buses, 1);
    
    for i = 1:n_buses
        % Calculate injected real and reactive power at each bus
        P_calc = 0; Q_calc = 0;
        
        for j = 1:n_buses
            if i ~= j
                P_calc = P_calc + abs(V(i)) * abs(V(j)) * (cos(theta(i) - theta(j)) * real(Ybus(i,j)) + sin(theta(i) - theta(j)) * imag(Ybus(i,j)));
                Q_calc = Q_calc + abs(V(i)) * abs(V(j)) * (cos(theta(i) - theta(j)) * imag(Ybus(i,j)) - sin(theta(i) - theta(j)) * real(Ybus(i,j)));
            end
        end
        
        % Power mismatch for bus i (P_mismatch, Q_mismatch)
        if bus_data(i, 2) == 3  % PQ bus
            P_mismatch(i) = P_calc - bus_data(i, 7);  % PL at PQ bus
            Q_mismatch(i) = Q_calc - bus_data(i, 8);  % QL at PQ bus
        elseif bus_data(i, 2) == 2  % PV bus
            P_mismatch(i) = P_calc - bus_data(i, 5);  % PG at PV bus
            Q_mismatch(i) = Q_calc - bus_data(i, 6);  % QG at PV bus
        end
    end
    
    % Check for convergence
    if max(abs([P_mismatch; Q_mismatch])) < tolerance
        converged = true;
        break;
    end
    
    % Step 5: Form Jacobian Matrix (J) for Newton-Raphson
    J = zeros(2*(n_buses - 1));  % Jacobian matrix size for voltage correction
    
    % Calculate partial derivatives for Jacobian (details omitted for brevity)
    % This part will involve calculating how voltage changes (dV, dθ) affect power mismatch
    
    % Solve for voltage corrections (dV, dθ)
    delta = J \ [P_mismatch(2:end); Q_mismatch(2:end)];
    
    % Update voltage magnitudes and angles
    for i = 2:n_buses
        theta(i) = theta(i) + delta(i - 1);               % Update angle (delta_theta)
        V(i) = V(i) + delta(n_buses - 1 + (i - 1));       % Update magnitude (delta_V)
    end
end

% Check for convergence
if converged
    disp(['Converged after ' num2str(iter) ' iterations']);
else
    disp('Did not converge within max iterations');
end

% Step 6: Calculate real and reactive power at each bus
for i = 1:n_buses
    P_calculated(i) = real(V(i)) * real(Ybus(i,i));
    Q_calculated(i) = imag(V(i)) * imag(Ybus(i,i));
end

% Step 7: Plot voltage profiles
figure;
subplot(2, 1, 1);
plot(1:n_buses, abs(V), '-o');
title('Voltage Magnitudes at Each Bus');
xlabel('Bus Number');
ylabel('Voltage (pu)');

subplot(2, 1, 2);
plot(1:n_buses, rad2deg(theta), '-o');
title('Voltage Angles at Each Bus');
xlabel('Bus Number');
ylabel('Voltage Angle (degrees)');
