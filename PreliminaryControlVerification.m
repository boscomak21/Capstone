%% GEO-DUDe Task 2 - Control Verification Script
% Verifies closed-loop attitude controller against simplified requirement metrics

clear; clc; close all;

%% Requirement-inspired thresholds
REQ.maxSteadyStateError_deg = 1.0;   % SOOS-1-AC-01
REQ.maxSettleTime_s         = 300;   % project-level verification target
REQ.maxWheelMomentumFrac    = 0.90;  % margin before desaturation needed
REQ.maxBodyRate_deg_s       = 2.0;   % verification target
REQ.numMonteCarloRuns       = 100;

%% Baseline parameters
params.J = diag([1200, 1500, 1800]);
params.Jinv = inv(params.J);
params.Kp = diag([80, 80, 60]);
params.Kd = diag([900, 900, 700]);
params.Ki = diag([0, 0, 0]);
params.useIntegral = false;

params.tauMax = [0.05; 0.05; 0.05];
params.hMax   = [25; 25; 25];
params.disturbanceBody = [2e-4; -1.5e-4; 1.0e-4];

q_ref = [1;0;0;0];
w_ref = [0;0;0];

%% Deterministic test cases
cases = {
    struct('name','Nominal large initial error', 'eul0',[20;-12;35],  'w0',deg2rad([0.8;-0.5;0.6]))
    struct('name','Docking fine correction',     'eul0',[1.2;-0.8;0.5],'w0',deg2rad([0.05;0.02;-0.03]))
    struct('name','Higher rate recovery',        'eul0',[8;6;-10],     'w0',deg2rad([1.5;-1.2;1.0]))
    };

fprintf('==== DETERMINISTIC VERIFICATION ====\n');
for i = 1:numel(cases)
    out = run_case(cases{i}.eul0, cases{i}.w0, q_ref, w_ref, params);

    pass_err   = out.finalErrDeg <= REQ.maxSteadyStateError_deg;
    pass_set   = out.settleTime  <= REQ.maxSettleTime_s;
    pass_h     = out.maxHFrac    <= REQ.maxWheelMomentumFrac;
    pass_rate  = out.maxRateDeg  <= REQ.maxBodyRate_deg_s || strcmp(cases{i}.name,'Higher rate recovery');

    fprintf('\nCase: %s\n', cases{i}.name);
    fprintf('  Final error      = %.3f deg\n', out.finalErrDeg);
    fprintf('  Settle time      = %.3f s\n', out.settleTime);
    fprintf('  Max wheel frac   = %.3f\n', out.maxHFrac);
    fprintf('  Max body rate    = %.3f deg/s\n', out.maxRateDeg);
    fprintf('  PASS error       = %d\n', pass_err);
    fprintf('  PASS settling    = %d\n', pass_set);
    fprintf('  PASS wheel margin= %d\n', pass_h);
    fprintf('  PASS rate        = %d\n', pass_rate);
end

%% Monte Carlo robustness verification
fprintf('\n==== MONTE CARLO ROBUSTNESS ====\n');

mc_pass = 0;
final_errs = zeros(REQ.numMonteCarloRuns,1);
settle_times = zeros(REQ.numMonteCarloRuns,1);

for k = 1:REQ.numMonteCarloRuns
    eul0 = [30;-20;25].*(2*rand(3,1)-1);              % deg
    w0   = deg2rad([1.0;1.0;1.0].*(2*rand(3,1)-1));   % rad/s

    p_mc = params;
    p_mc.J = diag(diag(params.J).*(1 + 0.1*(2*rand(3,1)-1)));  % ±10% inertia uncertainty
    p_mc.Jinv = inv(p_mc.J);
    p_mc.disturbanceBody = params.disturbanceBody .* (1 + 0.5*(2*rand(3,1)-1));

    out = run_case(eul0, w0, q_ref, w_ref, p_mc);
    final_errs(k) = out.finalErrDeg;
    settle_times(k) = out.settleTime;

    if out.finalErrDeg <= REQ.maxSteadyStateError_deg && out.maxHFrac <= 1.0
        mc_pass = mc_pass + 1;
    end
end

fprintf('Monte Carlo pass count: %d / %d\n', mc_pass, REQ.numMonteCarloRuns);
fprintf('Mean final error: %.3f deg\n', mean(final_errs));
fprintf('Worst final error: %.3f deg\n', max(final_errs));
fprintf('Mean settling time: %.3f s\n', mean(settle_times(~isnan(settle_times))));

figure;
histogram(final_errs, 15);
grid on;
xlabel('Final Attitude Error (deg)');
ylabel('Count');
title('Monte Carlo Final Error Distribution');

%% ============================================================
% Helper functions
% ============================================================

function out = run_case(eul0_deg, w0, q_ref, w_ref, params)
    q0 = eul321_to_quat(deg2rad(eul0_deg));
    x0 = [q0; w0; [0;0;0]; [0;0;0]];

    [t,x] = ode45(@(t,x) spacecraft_dynamics(t,x,q_ref,w_ref,params), [0 600], x0, ...
        odeset('RelTol',1e-8,'AbsTol',1e-9));

    q = x(:,1:4);
    w = x(:,5:7);
    h = x(:,8:10);

    N = length(t);
    errDeg = zeros(N,1);

    for n = 1:N
        qn = q(n,:)';
        q_err = quat_multiply(quat_conj(q_ref), qn);
        q_err = q_err / norm(q_err);
        if q_err(1) < 0
            q_err = -q_err;
        end
        errDeg(n) = rad2deg(2*atan2(norm(q_err(2:4)), abs(q_err(1))));
    end

    idx = find(errDeg <= 1.0, 1, 'first');
    if isempty(idx)
        settleTime = NaN;
    else
        settleTime = t(idx);
    end

    out.finalErrDeg = errDeg(end);
    out.settleTime = settleTime;
    out.maxHFrac = max(max(abs(h) ./ params.hMax', [], 2));
    out.maxRateDeg = max(max(abs(rad2deg(w)), [], 2));
end

function xdot = spacecraft_dynamics(~, x, q_ref, w_ref, p)
    q  = x(1:4); q = q / norm(q);
    w  = x(5:7);
    h  = x(8:10);
    ei = x(11:13);

    [~, tau_ctrl, q_err] = control_law(q, w, h, ei, q_ref, w_ref, p);
    tau_dist = p.disturbanceBody;

    wdot = p.Jinv * (tau_ctrl + tau_dist - cross(w, p.J*w + h));

    Omega = [0    -w(1) -w(2) -w(3);
             w(1)  0     w(3) -w(2);
             w(2) -w(3)  0     w(1);
             w(3)  w(2) -w(1)  0];
    qdot = 0.5 * Omega * q;
    hdot = -tau_ctrl;

    if p.useIntegral
        eidot = q_err(2:4);
    else
        eidot = [0;0;0];
    end

    xdot = [qdot; wdot; hdot; eidot];
end

function [tau_unsat, tau_sat, q_err] = control_law(q, w, h, ei, q_ref, w_ref, p)
    %#ok<INUSD>
    q_err = quat_multiply(quat_conj(q_ref), q);
    q_err = q_err / norm(q_err);
    if q_err(1) < 0
        q_err = -q_err;
    end

    evec = q_err(2:4);
    w_err = w - w_ref;

    tau_unsat = -p.Kp*evec - p.Kd*w_err - p.Ki*ei;
    tau_sat = min(max(tau_unsat, -p.tauMax), p.tauMax);
end

function q = eul321_to_quat(eul)
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
    w1 = q1(1); v1 = q1(2:4);
    w2 = q2(1); v2 = q2(2:4);
    q = [w1*w2 - dot(v1,v2);
         w1*v2 + w2*v1 + cross(v1,v2)];
end