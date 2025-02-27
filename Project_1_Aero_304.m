%% Project_1_Aero_304.m
% Circular Restricted Three-Body Problem (CR3BP) Analysis
% Yash Rotkar, Tricarico Jr. Robert Roland, Niko Angelakos
% Date: Feb 24th 2025

clc; 
clear;
close all; % Closes all figures if any

% Define constants
G = 1; % Gravitational constant (normalized)

%{
%% --- Part (a): Equations of Motion (Symbolic Verification - Optional) ---
% (This section is for symbolic verification and requires the Symbolic Math Toolbox)
% Commented out to prevent errors if the toolbox is not installed.
% syms x y xdot ydot mu real
% syms p1 p2 omega real
%
% % Kinetic Energy (T)
% % Potential Energy (V)
% % Lagrangian (L = T - V)
%
% % Derive equations of motion using Euler-Lagrange equations
%}

%% --- Part (b): Equilibrium Points (L4 and L5) ---
disp('--- Part (b): Equilibrium Points (L4 and L5) ---');
mu_values = [3.0039e-6, 1.2151e-2, 2.366e-4]; % Mass ratios for Sun-Earth, Earth-Moon, Saturn-Titan
system = {'Sun-Earth', 'Earth-Moon', 'Saturn-Titan'};

for i = 1:length(mu_values)
    mu = mu_values(i);

    % L4 coordinates (x, y)
    x_L4 = 0.5 - mu;
    y_L4 = sqrt(3)/2;

    % L5 coordinates (x, y)
    x_L5 = 0.5 - mu;
    y_L5 = -sqrt(3)/2;

    fprintf('System: %s\n', system{i});
    fprintf('  mu = %.6f\n', mu);
    fprintf('  L4: x = %.6f, y = %.6f\n', x_L4, y_L4);
    fprintf('  L5: x = %.6f, y = %.6f\n', x_L5, y_L5);
    fprintf('\n');
end

%% Part (c): Stability of Lagrange Points
disp('--- Part (c): Stability of Lagrange Points ---');
fprintf('\n');

% Define systems and mass ratios
syst = {'Sun-Earth', 'Earth-Moon', 'Saturn-Titan'};
mu_values = [3.0039e-6, 1.2151e-2, 2.366e-4]; % Given mass ratios

% Define second derivatives of U at Lagrange points
Uxx_val = [-2.99953, 8.862087, -1.0000056, 0.75, 0.75;  % Sun-Earth
           -2.84002, 7.38119, -1.0214, 0.75, 0.75;  % Earth-Moon
           -2.99224, 8.82602, -1.0004, 0.75, 0.75]; % Saturn-Titan

Uyy_val = [0, 0, 0, -1.25, -1.25;  % Sun-Earth
           0, 0, 0, -1.25, -1.25;  % Earth-Moon
           0, 0, 0, -1.25, -1.25]; % Saturn-Titan

Uxy_val = [0, 0, 0, 1.29903, -1.29903; % Sun-Earth
           0, 0, 0, 1.26747, -1.26747; % Earth-Moon
           0, 0, 0, 1.29838, -1.29838]; % Saturn-Titan

% Iterate over each system
for system = 1:length(syst)
    fprintf('\n     System: %s    \n', syst{system});

    % Iterates each Lagrange point
    for i = 1:5
        Uxx = Uxx_val(system, i);
        Uyy = Uyy_val(system, i);
        Uxy = Uxy_val(system, i);

        % Constructa the Jacobian matrix
        Y = [ 0   0   1  0;
              0   0   0  1;
             Uxx Uxy  0  2;
             Uxy Uyy -2  0];

        % Computes the eigenvalues
        eigenval = eig(Y);

        % Determines the stability Lagrange Point
        stable = all(real(eigenval) < 0.00001);
        stability_status = "Stable";
        fprintf('\n');
        
        if ~stable
            stability_status = "Unstable";
        end

        % Displays results
        fprintf('%s: Eigenvalues = [%.4f, %.4f, %.4f, %.4f]\n', ['L' num2str(i)], real(eigenval(1)), real(eigenval(2)), real(eigenval(3)), real(eigenval(4)));
        fprintf(' Stability: %s\n', stability_status);
    end
end


%{
 --- Part (c): Stability of Lagrange Points ---
disp('--- Part (c): Stability of Lagrange Points ---');
mu_values = [3.0039e-6, 1.2151e-2, 2.366e-4]; % Mass ratios
systems = {'Sun-Earth', 'Earth-Moon', 'Saturn-Titan'};

for i = 1:length(mu_values)
    mu = mu_values(i);
    fprintf('System: %s\n', systems{i});
    fprintf('  mu = %.6f\n', mu);
    fprintf('\n');

    % Iterate through Lagrange points L1 to L5 (approximate locations)
    if i == 1 %Sun-Earth
        L1 = 0.9900261;
        L2 = 1.0100345;
        L3 = -1.0000012;
    elseif i == 2 %Earth-Moon
        L1 = 0.836915;
        L2 = 1.15568;
        L3 = -1.00506;
    else %Saturn-Titan
        L1 = 0.9575;
        L2 = 1.0425;
        L3 = -1.0001;
    end

    Lagrange_points = [L1, L2, L3, 0.5-mu, 0.5-mu];
    for j = 1:5
        if j <=3
            x0 = Lagrange_points(j);
            y0 = 0; %L1,L2,L3
        elseif j == 4
            x0 = Lagrange_points(j);
            y0 = sqrt(3)/2;
        else
            x0 = Lagrange_points(j);
            y0 = -sqrt(3)/2;
        end

        % Compute Uxx, Uyy, Uxy numerically (perturbation method)
        delta = 1e-6; % Small perturbation
        Uxx = (U_potential(x0+delta, y0, mu) - 2*U_potential(x0, y0, mu) + U_potential(x0-delta, y0, mu)) / delta^2;
        Uyy = (U_potential(x0, y0+delta, mu) - 2*U_potential(x0, y0, mu) + U_potential(x0, y0-delta, mu)) / delta^2;
        Uxy = (U_potential(x0+delta, y0+delta, mu) - U_potential(x0+delta, y0-delta, mu) - U_potential(x0-delta, y0+delta, mu) + U_potential(x0-delta, y0-delta, mu)) / (4*delta^2);

        % Construct the A matrix
        A = [0 0 1 0;
             0 0 0 1;
             Uxx Uxy 0 2;
             Uxy Uyy -2 0];

        % Compute eigenvalues
        eigenvalues = eig(A);

        fprintf('  Lagrange point L%d: x = %.6f, y = %.6f\n', j, x0, y0);
        fprintf('    Eigenvalues: %.6f %.6f %.6f %.6f\n', eigenvalues(1), eigenvalues(2), eigenvalues(3), eigenvalues(4));

        % Stability assessment (simple)
        if all(real(eigenvalues) == 0) && all(imag(eigenvalues) ~= 0)
            fprintf('    ==> LINEARLY STABLE (center)\n');
            fprintf('\n');
        else
            fprintf('    ==> UNSTABLE (saddle or spiral)\n');
            fprintf('\n');
        end
    end
    fprintf('\n');
end
%}


%{
%% --- Part (d): Simulating Orbits in the Earth-Moon System around L2 point ---
disp('--- Part (d): Simulating Orbits around L2 ---');
load('EM_L2-304P1.mat'); % Loads x0, T, mu, perturbation

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Numerical integration for the nominal orbit
[t_nominal, xy_nominal] = ode45(@(t,xy) cr3bp_eom(t,xy,mu), [0, T], x0, options);

% Numerical integration for the perturbed orbit
x0_perturbed = x0 + perturbation;
[t_perturbed, xy_perturbed] = ode45(@(t,xy) cr3bp_eom(t,xy,mu), [0, T], x0_perturbed, options);

% Plot nominal orbit in rotating frame
figure;
plot(xy_nominal(:,1), xy_nominal(:,2));
xlabel('x');
ylabel('y');
title('L2 Nominal Orbit in Rotating Frame');
axis equal;
grid on;

% Convert to inertial frame
X_N = zeros(size(xy_nominal, 1), 1);
Y_N = zeros(size(xy_nominal, 1), 1);

for i = 1:length(t_nominal)
    theta = t_nominal(i); % Assuming omega = 1
    X_N(i) = xy_nominal(i,1)*cos(theta) - xy_nominal(i,2)*sin(theta);
    Y_N(i) = xy_nominal(i,1)*sin(theta) + xy_nominal(i,2)*cos(theta);
end

figure;
plot(X_N, Y_N);
xlabel('X');
ylabel('Y');
title('L2 Nominal Orbit in Inertial Frame');
axis equal;
grid on;

% Compute departure motion
delta_xy = xy_perturbed - xy_nominal;

% Plot departure motion
figure;
subplot(2,1,1);
plot(t_perturbed, delta_xy(:,1:2));
xlabel('Time');
ylabel('delta x, delta y');
legend('delta x', 'delta y');
title('L2 Departure Motion (Position)');
grid on;

subplot(2,1,2);
plot(t_perturbed, delta_xy(:,3:4));
xlabel('Time');
ylabel('delta xdot, delta ydot');
legend('delta xdot', 'delta ydot');
title('L2 Departure Motion (Velocity)');
grid on;

% Linearized equations of motion
% Compute Uxx, Uyy, Uxy numerically (perturbation method)

delta = 1e-6; % Small perturbation

Uxx = (U_potential(x0_perturbed(1)+delta, x0_perturbed(2), mu) - 2*U_potential(x0_perturbed(1), x0_perturbed(2), mu) + U_potential(x0_perturbed(1)-delta, x0_perturbed(2), mu)) / delta^2;
Uyy = (U_potential(x0_perturbed(1), x0_perturbed(2)+delta, mu) - 2*U_potential(x0_perturbed(1), x0_perturbed(2), mu) + U_potential(x0_perturbed(1), x0_perturbed(2)-delta, mu)) / delta^2;
Uxy = (U_potential(x0_perturbed(1)+delta, x0_perturbed(2)+delta, mu) - U_potential(x0_perturbed(1)+delta, x0_perturbed(2)-delta, mu) - U_potential(x0_perturbed(1)-delta, x0_perturbed(2)+delta, mu) + U_potential(x0_perturbed(1)-delta, x0_perturbed(2)-delta, mu)) / (4*delta^2);

A = [0 0 1 0;
     0 0 0 1;
     Uxx Uxy 0 2;
     Uxy Uyy -2 0];

% Integrate the linearized equation of motion
delta_x0 = perturbation;
[t_linear, delta_xy_linear] = ode45(@(t,delta_xy) A*delta_xy, [0, T], delta_x0, options);

% Compare with nonlinear departure motion
figure;
subplot(2,1,1);
plot(t_perturbed, delta_xy(:,1:2), t_linear, delta_xy_linear(:,1:2));
xlabel('Time');
ylabel('delta x, delta y');
legend('Nonlinear delta x', 'Nonlinear delta y', 'Linear delta x', 'Linear delta y');
title('L2 Comparison of Departure Motion (Position)');
grid on;

subplot(2,1,2);
plot(t_perturbed, delta_xy(:,3:4), t_linear, delta_xy_linear(:,3:4));
xlabel('Time');
ylabel('delta xdot, delta ydot');
legend('Nonlinear delta xdot', 'Nonlinear delta ydot', 'Linear delta xdot', 'Linear delta ydot');
title('L2 Comparison of Departure Motion (Velocity)');
grid on;

% --- Part (e): Simulating Orbits in the Earth-Moon System around L4 point ---
disp('--- Part (e): Simulating Orbits around L4 ---');
load('EM_L4-304P1.mat'); % Loads x0, T, mu, perturbation

% Numerical integration for the nominal orbit
[t_nominal, xy_nominal] = ode45(@(t,xy) cr3bp_eom(t,xy,mu), [0, T], x0, options);

% Numerical integration for the perturbed orbit
x0_perturbed = x0 + perturbation;
[t_perturbed, xy_perturbed] = ode45(@(t,xy) cr3bp_eom(t,xy,mu), [0, T], x0_perturbed, options);

% Plot nominal orbit in rotating frame
figure;
plot(xy_nominal(:,1), xy_nominal(:,2));
xlabel('x');
ylabel('y');
title('L4 Nominal Orbit in Rotating Frame');
axis equal;
grid on;

% Convert to inertial frame
X_N = zeros(size(xy_nominal, 1), 1);
Y_N = zeros(size(xy_nominal, 1), 1);

for i = 1:length(t_nominal)
    theta = t_nominal(i); % Assuming omega = 1
    X_N(i) = xy_nominal(i,1)*cos(theta) - xy_nominal(i,2)*sin(theta);
    Y_N(i) = xy_nominal(i,1)*sin(theta) + xy_nominal(i,2)*cos(theta);
end

figure;
plot(X_N, Y_N);
xlabel('X');
ylabel('Y');
title('L4 Nominal Orbit in Inertial Frame');
axis equal;
grid on;

% Compute departure motion
delta_xy = xy_perturbed - xy_nominal;

% Plot departure motion
figure;
subplot(2,1,1);
plot(t_perturbed, delta_xy(:,1:2));
xlabel('Time');
ylabel('delta x, delta y');
legend('delta x', 'delta y');
title('L4 Departure Motion (Position)');
grid on;

subplot(2,1,2);
plot(t_perturbed, delta_xy(:,3:4));
xlabel('Time');
ylabel('delta xdot, delta ydot');
legend('delta xdot', 'delta ydot');
title('L4 Departure Motion (Velocity)');
grid on;

% Linearized equations of motion
% Compute Uxx, Uyy, Uxy numerically (perturbation method)

delta = 1e-6; % Small perturbation
Uxx = (U_potential(x0_perturbed(1)+delta, x0_perturbed(2), mu) - 2*U_potential(x0_perturbed(1), x0_perturbed(2), mu) + U_potential(x0_perturbed(1)-delta, x0_perturbed(2), mu)) / delta^2;
Uyy = (U_potential(x0_perturbed(1), x0_perturbed(2)+delta, mu) - 2*U_potential(x0_perturbed(1), x0_perturbed(2), mu) + U_potential(x0_perturbed(1), x0_perturbed(2)-delta, mu)) / delta^2;
Uxy = (U_potential(x0_perturbed(1)+delta, x0_perturbed(2)+delta, mu) - U_potential(x0_perturbed(1)+delta, x0_perturbed(2)-delta, mu) - U_potential(x0_perturbed(1)-delta, x0_perturbed(2)+delta, mu) + U_potential(x0_perturbed(1)-delta, x0_perturbed(2)-delta, mu)) / (4*delta^2);

A = [0 0 1 0;
     0 0 0 1;
     Uxx Uxy 0 2;
     Uxy Uyy -2 0];

% Integrate the linearized equation of motion
delta_x0 = perturbation;
[t_linear, delta_xy_linear] = ode45(@(t,delta_xy) A*delta_xy, [0, T], delta_x0, options);

% Compare with nonlinear departure motion
figure;
subplot(2,1,1);
plot(t_perturbed, delta_xy(:,1:2), t_linear, delta_xy_linear(:,1:2));
xlabel('Time');
ylabel('delta x, delta y');
legend('Nonlinear delta x', 'Nonlinear delta y', 'Linear delta x', 'Linear delta y');
title('L4 Comparison of Departure Motion (Position)');
grid on;

subplot(2,1,2);
plot(t_perturbed, delta_xy(:,3:4), t_linear, delta_xy_linear(:,3:4));
xlabel('Time');
ylabel('delta xdot', 'delta ydot');
legend('Nonlinear delta xdot', 'Nonlinear delta ydot', 'Linear delta xdot', 'Linear delta ydot');
title('L4 Comparison of Departure Motion (Velocity)');
grid on;

disp('Simulation complete.');

% --- Function Definitions ---
% CR3BP Equations of Motion (normalized)

 function dxydt = cr3bp_eom(t, xy, mu)
    x = xy(1);
    y = xy(2);
    xdot = xy(3);
    ydot = xy(4);

    p1 = sqrt((x + mu)^2 + y^2);
    p2 = sqrt((x - 1 + mu)^2 + y^2);

    xddot = x + 2*ydot - (1-mu)*(x+mu)/p1^3 - mu*(x-1+mu)/p2^3;
    yddot = y - 2*xdot - (1-mu)*y/p1^3 - mu*y/p2^3;

    dxydt = [xdot; ydot; xddot; yddot];
end

% Effective Potential
function U = U_potential(x, y, mu)
    p1 = sqrt((x + mu)^2 + y^2);
    p2 = sqrt((x - 1 + mu)^2 + y^2);
    U = 0.5*(x^2 + y^2) + (1-mu)/p1 + mu/p2;
end
%}

