%% GEO-DUDe Task 2 - 3-Axis Quaternion Attitude Control
% Quaternion-based PD/PID attitude control with reaction wheel torque limits
% Suitable as an early WP2107 control-law baseline for closed-loop simulation

clear; clc; close all;

%% ---------------------------
%  Spacecraft / control params
% ----------------------------
params.J = diag([1200, 1500, 1800]);   % kg*m^2, example servicer principal inertias
params.Jinv = inv(params.J);

% Control gains (tune these)
params.Kp = diag([80, 80, 60]);        % proportional gain on quaternion vector error
params.Kd = diag([900, 900, 700]);     % derivative gain on body rate
params.Ki = diag([0, 0, 0]);           % start with zero; enable later if needed

% Wheel / actuator limits
params.tauMax = [0.05; 0.05; 0.05];    % N*m per axis
params.hMax   = [25; 25; 25];          % N*m*s wheel momentum capacity
params.useIntegral = false;

% Simple disturbance model (gravity gradient / plume / SRP lump)
params.disturbanceBody = [2e-4; -1.5e-4; 1.0e-4];  % N*m constant bias

% Reference attitude and rate
q_ref = [1; 0; 0; 0];   % desired quaternion [qw qx qy qz]'
w_ref = [0; 0; 0];      % rad/s

%% ---------------------------
%  Initial conditions
% ----------------------------
eul0_deg = [20; -12; 35];                 % initial yaw/pitch/roll style guess
q0 = eul321_to_quat(deg2rad(eul0_deg));   % initial attitude
w0 = deg2rad([0.8; -0.5; 0.6]);           % rad/s initial body rates
h0 = [0; 0; 0];                           % wheel momentum
ei0 = [0; 0; 0];                          % integral of attitude error

x0 = [q0; w0; h0; ei0];

%% ---------------------------
%  Simulation settings
% ----------------------------
tspan = [0 600];   % s
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);

[t, x] = ode45(@(t,x) spacecraft_dynamics(t, x, q_ref, w_ref, params), tspan, x0, opts);

%% ---------------------------
%  Post-process
% ----------------------------
N = length(t);
q = x(:,1:4);
w = x(:,5:7);
h = x(:,8:10);
ei = x(:,11:13); %#ok<NASGU>

theta_err_deg = zeros(N,1);
tau_cmd = zeros(N,3);
tau_applied = zeros(N,3);

for k = 1:N
    qk = q(k,:)';
    wk = w(k,:)';
    hk = h(k,:)';
    eik = x(k,11:13)';

    [tau_u, tau_sat, q_err] = control_law(qk, wk, hk, eik, q_ref, w_ref, params);
    tau_cmd(k,:) = tau_u';
    tau_applied(k,:) = tau_sat';
    theta_err_deg(k) = rad2deg(2*atan2(norm(q_err(2:4)), abs(q_err(1))));
end

settle_idx = find(theta_err_deg <= 1.0, 1, 'first');
if isempty(settle_idx)
    settle_time = NaN;
else
    settle_time = t(settle_idx);
end

peak_err = max(theta_err_deg);
final_err = theta_err_deg(end);

fprintf('Peak attitude error: %.3f deg\n', peak_err);
fprintf('Final attitude error: %.3f deg\n', final_err);
fprintf('Time to enter 1 deg band: %.3f s\n', settle_time);

%% ---------------------------
%  Plots
% ----------------------------
figure;
plot(t, theta_err_deg, 'LineWidth', 1.5); hold on;
yline(1.0, '--r', '1 deg requirement');
grid on;
xlabel('Time (s)');
ylabel('Attitude Error (deg)');
title('GEO-DUDe Attitude Error');

figure;
plot(t, rad2deg(w), 'LineWidth', 1.3);
grid on;
xlabel('Time (s)');
ylabel('Body Rate (deg/s)');
legend('\omega_x','\omega_y','\omega_z','Location','best');
title('Body Angular Rates');

figure;
plot(t, tau_applied, 'LineWidth', 1.3); hold on;
yline(params.tauMax(1), '--k');
yline(-params.tauMax(1), '--k');
grid on;
xlabel('Time (s)');
ylabel('Applied Torque (N·m)');
legend('\tau_x','\tau_y','\tau_z','Location','best');
title('Applied Control Torque');

figure;
plot(t, h, 'LineWidth', 1.3); hold on;
yline(params.hMax(1), '--r');
yline(-params.hMax(1), '--r');
grid on;
xlabel('Time (s)');
ylabel('Wheel Momentum (N·m·s)');
legend('h_x','h_y','h_z','Location','best');
title('Reaction Wheel Momentum');

%% ============================================================
% Local functions
% ============================================================

function xdot = spacecraft_dynamics(~, x, q_ref, w_ref, p)
    q  = x(1:4);
    w  = x(5:7);
    h  = x(8:10);
    ei = x(11:13);

    q = q / norm(q);

    [~, tau_ctrl, q_err] = control_law(q, w, h, ei, q_ref, w_ref, p);

    % External disturbance torque
    tau_dist = p.disturbanceBody;

    % Rigid-body rotational dynamics: J*w_dot + w×(Jw + h) = tau_ctrl + tau_dist
    wdot = p.Jinv * (tau_ctrl + tau_dist - cross(w, p.J*w + h));

    % Quaternion kinematics
    Omega = [0    -w(1) -w(2) -w(3);
             w(1)  0     w(3) -w(2);
             w(2) -w(3)  0     w(1);
             w(3)  w(2) -w(1)  0];
    qdot = 0.5 * Omega * q;

    % Wheel momentum accumulation (simple model)
    hdot = -tau_ctrl;

    % Integral of vector attitude error
    if p.useIntegral
        evec = q_err(2:4);
        eidot = evec;
    else
        eidot = [0;0;0];
    end

    xdot = [qdot; wdot; hdot; eidot];
end

function [tau_unsat, tau_sat, q_err] = control_law(q, w, h, ei, q_ref, w_ref, p)
    %#ok<INUSD> h currently included only for monitoring / future desat logic

    % Quaternion error from current to reference
    q_err = quat_multiply(quat_conj(q_ref), q);
    q_err = q_err / norm(q_err);

    % Ensure shortest rotation
    if q_err(1) < 0
        q_err = -q_err;
    end

    evec = q_err(2:4);
    w_err = w - w_ref;

    tau_unsat = -p.Kp*evec - p.Kd*w_err - p.Ki*ei;

    tau_sat = min(max(tau_unsat, -p.tauMax), p.tauMax);
end

function q = eul321_to_quat(eul)
    % yaw-pitch-roll: psi theta phi
    psi = eul(1); th = eul(2); phi = eul(3);
    c1 = cos(psi/2); s1 = sin(psi/2);
    c2 = cos(th/2);  s2 = sin(th/2);
    c3 = cos(phi/2); s3 = sin(phi/2);

    q = [ c1*c2*c3 + s1*s2*s3;
          c1*c2*s3 - s1*s2*c3;
          c1*s2*c3 + s1*c2*s3;
          s1*c2*c3 - c1*s2*s3 ];
    q = q / norm(q);
end

function qc = quat_conj(q)
    qc = [q(1); -q(2:4)];
end

function q = quat_multiply(q1, q2)
    % Hamilton product, scalar-first
    w1 = q1(1); v1 = q1(2:4);
    w2 = q2(1); v2 = q2(2:4);

    q = [w1*w2 - dot(v1,v2);
         w1*v2 + w2*v1 + cross(v1,v2)];
end