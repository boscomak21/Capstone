%% Cube face pointing (single-axis attitude) with PID control
% Graph angular position, velocity, and acceleration as the cube rotates to a setpoint.
% Units: rad, rad/s, rad/s^2

clear; clc; close all;

%% --- Parameters you can tweak ---
I = 0.02;                 % kg*m^2  (effective moment of inertia about the axis)
tauMax = 0.05;            % N*m     (actuator torque limit)
zeta = 0.9;               % desired damping ratio (0.7–1.2 typical)
wn   = 2.5;               % desired natural frequency rad/s (bigger = faster)
usePID = true;            % set false to use PD only

% Setpoint (target angle for the cube face about this axis)
theta_sp = deg2rad(45);   % rad

% Arbitrary initial conditions (position, velocity). Accel is implied by dynamics.
theta0 = deg2rad(-110);   % rad
omega0 = deg2rad(40);     % rad/s

% Optional: a constant disturbance torque (e.g., gravity gradient / SRP lumped)
tauDist = 0.0;            % N*m

% Simulation time
tEnd = 10;                % s

%% --- Controller gains (designed from 2nd-order target behavior) ---
% For rigid-body single-axis: I*theta_ddot = tau
% PD gives closed-loop approx: theta_ddot + (Kd/I)*theta_dot + (Kp/I)*theta = (Kp/I)*theta_sp
Kp = I * wn^2;
Kd = 2 * zeta * wn * I;

% Add integral action if requested (helps eliminate steady-state error under disturbance)
Ki = 0;
if usePID
    % A reasonable starting point for Ki (tune if it overshoots/oscillates):
    Ki = 0.35 * Kp * wn;  % heuristic
end

%% --- Simulate with ODE45 ---
% State x = [theta; omega; integral_error]
x0 = [theta0; omega0; 0];

opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
[t, x] = ode45(@(t,x) dynamics(t, x, I, Kp, Kd, Ki, tauMax, tauDist, theta_sp), [0 tEnd], x0, opts);

theta = x(:,1);
omega = x(:,2);

% Compute commanded torque and acceleration for plotting
tauCmd = zeros(size(t));
alpha  = zeros(size(t));
for i = 1:numel(t)
    [~, tauCmd(i), alpha(i)] = dynamics(t(i), x(i,:).', I, Kp, Kd, Ki, tauMax, tauDist, theta_sp);
end

%% --- Plot ---
figure;
subplot(3,1,1);
plot(t, rad2deg(theta), 'LineWidth', 1.5); hold on;
yline(rad2deg(theta_sp), '--', 'Setpoint');
grid on;
ylabel('\theta (deg)');
title('Cube Face Pointing: Angular Position / Velocity / Acceleration');

subplot(3,1,2);
plot(t, rad2deg(omega), 'LineWidth', 1.5);
grid on;
ylabel('\omega (deg/s)');

subplot(3,1,3);
plot(t, rad2deg(alpha), 'LineWidth', 1.5); hold on;
grid on;
ylabel('\alpha (deg/s^2)');
xlabel('Time (s)');

figure;
plot(t, tauCmd, 'LineWidth', 1.5); hold on;
yline(tauMax,'--','+Torque limit'); yline(-tauMax,'--','-Torque limit');
grid on;
xlabel('Time (s)'); ylabel('\tau (N·m)');
title('Commanded Control Torque');

%% --- Dynamics function (local) ---
function [xdot, tau, alpha] = dynamics(~, x, I, Kp, Kd, Ki, tauMax, tauDist, theta_sp)
    theta = x(1);
    omega = x(2);
    eint  = x(3);

    % Angle error (wrapped to [-pi, pi] so it takes the shortest way)
    e = wrapToPi(theta_sp - theta);

    % PID torque command
    tau_unsat = Kp*e + Kd*(0 - omega) + Ki*eint;

    % Saturate actuator torque
    tau = max(min(tau_unsat, tauMax), -tauMax);

    % Rigid-body rotational dynamics
    % I*theta_ddot = tau + disturbance
    alpha = (tau + tauDist) / I;

    % State derivatives
    xdot = [omega; alpha; e];
end

function a = wrapToPi(a)
    a = mod(a + pi, 2*pi) - pi;
end