function figure_total_blast_vs_lambda_per_parameter()
% One curve per parameter: effect of scaling lambda on AUC = int_0^T B(t) dt
% T = 90 days (3 months)
%
% Curves:
%   k1 (black)
%   initial response group (GREEN): N, gamma, k2, E0, k4   (each separate curve)
%   long term group (BLUE): delta, eps, k3                (each separate curve)
%   B0 (red)
%
% Lambda sampling:
%   - For all parameters except N: lambda = logspace(lamMinPow, lamMaxPow, nPts)
%   - For N only: lambda = N/N0 for integer N (within range)

clearvars; 
% close all; 
clc;

%% ---------------------- settings ----------------------
Tend = 90;
lamMin = 1e-2;
lamMax = 1e1;
nPts   = 31;                          % sampling for non-N parameters
lamVec = logspace(log10(lamMin), log10(lamMax), nPts);

odeOpts = odeset('RelTol',1e-7,'AbsTol',1e-10);

%% ---------------------- baseline model ----------------------
P0  = default_params();
IC0 = default_initial_conditions(P0);
N0  = P0.N;

E0_total0 = sum(IC0.E0);

%% ---------------------- lambda values for N (integer constraint) ----------------------
Nmin = max(2, ceil(lamMin * N0));
Nmax = max(Nmin, floor(lamMax * N0));
Nvals  = Nmin:Nmax;
lamVecN = Nvals / N0;

%% ---------------------- colors (match your figure intent) ----------------------
col_k1   = [0 0 0];          % black
col_init = [0.1 0.5 0.1];    % green-ish
col_long = [0 0.2 0.9];      % blue-ish
col_B0   = [0.9 0 0];        % red-ish

%% ---------------------- compute curves ----------------------
% Helper anonymous for AUC
AUC = @(P,IC) blast_auc(P, IC, Tend, odeOpts);

totalJobs = numel(lamVec)*9 + numel(lamVecN);
done = 0;
t0 = tic;

progressPrint = @(done) fprintf('\rProgress: %d/%d (%.1f%%) | Elapsed: %s', ...
    done, totalJobs, 100*done/totalJobs, prettyDuration(toc(t0)));

% k1 curve (black)
Y_k1 = zeros(size(lamVec));
for j = 1:numel(lamVec)
    lam = lamVec(j);
    P = P0; IC = IC0;
    P.k1 = P0.k1 * lam;
    Y_k1(j) = AUC(P,IC);
    done = done + 1;
    if mod(done,10)==0 || done==totalJobs
        progressPrint(done);
    end
end

% INITIAL RESPONSE: gamma (green)
Y_gamma = zeros(size(lamVec));
for j = 1:numel(lamVec)
    lam = lamVec(j);
    P = P0; IC = IC0;
    P.gamma = P0.gamma * lam;
    Y_gamma(j) = AUC(P,IC);
    if mod(done,10)==0 || done==totalJobs
        progressPrint(done);
    end
end

% Switching threshold: Bhalf (I suggest grouping with initial response, but your choice)
Y_Bhalf = zeros(size(lamVec));
for j = 1:numel(lamVec)
    lam = lamVec(j);
    P = P0; IC = IC0;
    P.Bhalf = P0.Bhalf * lam;
    Y_Bhalf(j) = AUC(P,IC);
    if mod(done,10)==0 || done==totalJobs
        progressPrint(done);
    end
end

% INITIAL RESPONSE: k2 (green)
Y_k2 = zeros(size(lamVec));
for j = 1:numel(lamVec)
    lam = lamVec(j);
    P = P0; IC = IC0;
    P.k2 = P0.k2 * lam;
    Y_k2(j) = AUC(P,IC);
    if mod(done,10)==0 || done==totalJobs
        progressPrint(done);
    end
end

% INITIAL RESPONSE: k4 (green)
Y_k4 = zeros(size(lamVec));
for j = 1:numel(lamVec)
    lam = lamVec(j);
    P = P0; IC = IC0;
    P.k4 = P0.k4 * lam;
    Y_k4(j) = AUC(P,IC);
    if mod(done,10)==0 || done==totalJobs
        progressPrint(done);
    end
end

% % INITIAL RESPONSE: E0 (green) -> scale total initial effector pool only
% Y_E0 = zeros(size(lamVec));
% for j = 1:numel(lamVec)
%     lam = lamVec(j);
%     P = P0;
%     IC = IC0;
%     E0_total = E0_total0 * lam;
%     IC.E0 = (E0_total / P.N) * ones(P.N,1);
%     Y_E0(j) = AUC(P,IC);
%     if mod(done,10)==0 || done==totalJobs
%         progressPrint(done);
%     end
% end

% INITIAL RESPONSE: CAR0 (green) -> scale total CAR T cells, keep partition fixed
Y_CAR0 = zeros(size(lamVec));
CAR0_base = (IC0.M0 / 0.468);  % recover baseline CAR0 from IC (consistent with your IC builder)

for j = 1:numel(lamVec)
    lam = lamVec(j);

    P  = P0;
    IC = IC0;

    CAR0 = CAR0_base * lam;

    % keep partition fixed (same as your default_initial_conditions)
    IC.M0 = 0.468 * CAR0;
    IC.E0 = (0.532 * CAR0 / P.N) * ones(P.N,1);

    Y_CAR0(j) = AUC(P, IC);
end

% INITIAL RESPONSE: N (green) -> N must be integer; keep total E0 constant and redistribute
Y_N = zeros(size(lamVecN));
for j = 1:numel(lamVecN)
    P = P0;
    P.N = Nvals(j);

    IC = default_initial_conditions(P);  % correct vector sizes
    % Keep same baseline totals for B0, A0, M0
    IC.B0 = IC0.B0;
    IC.A0 = IC0.A0;
    IC.M0 = IC0.M0;

    % Keep total effector pool constant, just redistributed across N bins
    IC.E0 = (E0_total0 / P.N) * ones(P.N,1);

    Y_N(j) = AUC(P,IC);
    if mod(done,10)==0 || done==totalJobs
        progressPrint(done);
    end
end

% LONG TERM: delta (blue)
Y_delta = zeros(size(lamVec));
for j = 1:numel(lamVec)
    lam = lamVec(j);
    P = P0; IC = IC0;
    P.delta = P0.delta * lam;
    Y_delta(j) = AUC(P,IC);
    if mod(done,10)==0 || done==totalJobs
        progressPrint(done);
    end
end

% LONG TERM: eps (blue)
Y_eps = zeros(size(lamVec));
for j = 1:numel(lamVec)
    lam = lamVec(j);
    P = P0; IC = IC0;
    P.eps = P0.eps * lam;
    Y_eps(j) = AUC(P,IC);
    if mod(done,10)==0 || done==totalJobs
        progressPrint(done);
    end
end

% LONG TERM: k3 (blue)
Y_k3 = zeros(size(lamVec));
for j = 1:numel(lamVec)
    lam = lamVec(j);
    P = P0; IC = IC0;
    P.k3 = P0.k3 * lam;
    Y_k3(j) = AUC(P,IC);
    if mod(done,10)==0 || done==totalJobs
        progressPrint(done);
    end
end

% B0 (red)
Y_B0 = zeros(size(lamVec));
for j = 1:numel(lamVec)
    lam = lamVec(j);
    P = P0; IC = IC0;
    IC.B0 = IC0.B0 * lam;
    Y_B0(j) = AUC(P,IC);
    if mod(done,10)==0 || done==totalJobs
        progressPrint(done);
    end
end

fprintf('\n');

%% ---------------------- plot (log-log) ----------------------
fig = figure('Color','w','Position',[120 120 900 560]);
ax = axes(fig); hold(ax,'on');
Y_k1    = max(Y_k1,    1e-12);
Y_gamma = max(Y_gamma, 1e-12);
Y_k2    = max(Y_k2,    1e-12);
Y_k4    = max(Y_k4,    1e-12);
% Y_E0    = max(Y_E0,    1e-12);
Y_CAR0 = max(Y_CAR0, 1e-12);
Y_N     = max(Y_N,     1e-12);
Y_Bhalf = max(Y_Bhalf, 1e-12);
Y_delta = max(Y_delta, 1e-12);
Y_eps   = max(Y_eps,   1e-12);
Y_k3    = max(Y_k3,    1e-12);
Y_B0    = max(Y_B0,    1e-12);
% k1 (black)
loglog(ax, lamVec, Y_k1, '-o', 'LineWidth', 2, 'Color', col_k1, 'MarkerSize', 5);

% initial response (green) - separate curves, same color
loglog(ax, lamVec,  Y_gamma, '-o', 'LineWidth', 2, 'Color', col_init, 'MarkerSize', 5);
loglog(ax, lamVec,  Y_k2,    '--', 'LineWidth', 2, 'Color', col_init, 'MarkerSize', 5);
loglog(ax, lamVec,  Y_k4,    '-x', 'LineWidth', 2, 'Color', col_init, 'MarkerSize', 5);
% loglog(ax, lamVec,  Y_E0,    '-.', 'LineWidth', 2, 'Color', col_init, 'MarkerSize', 5);
loglog(ax, lamVec, Y_CAR0, '-.', 'LineWidth', 2, 'Color', col_init, 'MarkerSize', 5);
loglog(ax, lamVecN, Y_N,     ':', 'LineWidth', 2, 'Color', col_init, 'MarkerSize', 5);

% long term (blue) - separate curves, same color
loglog(ax, lamVec, Y_Bhalf, '--', 'LineWidth', 2, 'Color', col_long, 'MarkerSize', 5);
loglog(ax, lamVec, Y_delta, '-o', 'LineWidth', 2, 'Color', col_long, 'MarkerSize', 5);
loglog(ax, lamVec, Y_eps,   '-x', 'LineWidth', 2, 'Color', col_long, 'MarkerSize', 5);
loglog(ax, lamVec, Y_k3,    '-.', 'LineWidth', 2, 'Color', col_long, 'MarkerSize', 5);

% B0 (red)
loglog(ax, lamVec, Y_B0, '-o', 'LineWidth', 2, 'Color', col_B0, 'MarkerSize', 5);

grid(ax,'on'); box(ax,'on');

xlabel(ax, '$\theta$, scaling parameter', 'Interpreter','latex', 'FontSize',22);
ylabel(ax, '$\int_0^T B(t)\,dt$', 'Interpreter','latex', 'FontSize',22);
title(ax, 'Total blast cells up to 3 months', 'Interpreter','latex', 'FontSize',25);
xlim(ax, [lamMin lamMax]);

% Legend: one entry per parameter curve


% lg = legend(ax, ...
%     {'k_1', ...
%      '\gamma', 'k_2', 'k_4', 'E_0 (total)', 'N (integer)', 'B_{1/2}', ...
%      '\delta', '\epsilon', 'k_3', ...
%      'B_0'}, ...
%     'Location','southeast');
% lg = legend(ax, ...
%     {'k_1', ...
%      '\gamma', 'k_2', 'k_4', 'CAR_0 (dose)', 'N (integer)', 'B_{1/2}', ...
%      '\delta', '\epsilon', 'k_3', ...
%      'B_0'}, ...
%     'Location','southeast');
% set(lg, 'Box','on', 'LineWidth', 2);
set(gca,'XScale','log','YScale','log');
set(ax,'XScale','log','YScale','log');

%% ---------------------- 3-Column Custom Legend (Ultra-Tight) ----------------------
% Create a wide, short legend box at the bottom
% Position: [left, bottom, width, height]
% lg_ax = axes('Position', [0.12, 0.11, 0.80, 0.18], 'Units', 'normalized');
lg_ax = axes('Position', [0.43, 0.12, 0.46, 0.25], 'Units', 'normalized'); 
hold(lg_ax, 'on');
set(lg_ax, 'XTick', [], 'YTick', [], 'Box', 'on', 'LineWidth', 1.2, 'Color', 'w');

% Layout Settings - TIGHTENED
% We reduce the gap between columns and the gap between text and lines
col_x = [0.01, 0.33, 0.70]; % Starting X for the 3 columns (closer together)
line_off = 0.1;            % MUCH TIGHTER: Distance from text start to line start
line_len = 0.08;            % Standard line swatch length
y_top = 0.88;
y_step = 0.155;              % Vertical spacing

% Helper to draw a header
draw_h = @(txt, col, xp, yp) text(xp, yp, txt, 'Color', col, ... ...
    'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'latex', 'Parent', lg_ax);

% Helper to draw a parameter row
    function draw_p(txt, col, style, mark, xp, yp)
        % 1. Place text
        text(xp + 0.01, yp, txt, 'Color', 'k', 'Interpreter', 'latex', ...
             'FontSize', 15, 'Parent', lg_ax);
        
        % 2. Define 3 points for the line [start, middle, end]
        lx = [xp + line_off, xp + line_off + line_len/2, xp + line_off + line_len];
        ly = [yp, yp, yp];
        
        % 3. Plot line with centered marker
        plot(lg_ax, lx, ly, 'Color', col, 'LineStyle', style, 'LineWidth', 2, ...
            'Marker', mark, 'MarkerIndices', 2, ... % Forces marker only on the middle point
            'MarkerSize', 5);
    end
% --- COLUMN 1: Initial Response ---
curr_y = y_top;
draw_h('Initial Response', col_init, col_x(1), curr_y);
curr_y = curr_y - y_step;
draw_p('$N$', col_init, ':', 'none', col_x(1), curr_y); curr_y = curr_y - y_step;
draw_p('$\gamma$', col_init, '-', 'o', col_x(1), curr_y); curr_y = curr_y - y_step;
draw_p('$k_2$', col_init, '--', 'none', col_x(1), curr_y); curr_y = curr_y - y_step;
draw_p('$k_4$', col_init, '-', 'x', col_x(1), curr_y); curr_y = curr_y - y_step;
draw_p('$E_0$', col_init, '-.', 'none', col_x(1), curr_y); 

% --- COLUMN 2: Long Term Response ---
curr_y = y_top;
draw_h('Long Term Response', col_long, col_x(2), curr_y);
curr_y = curr_y - y_step;
draw_p('$B_{1/2}$', col_long, '--', 'none', col_x(2), curr_y); curr_y = curr_y - y_step;
draw_p('$\delta$', col_long, '-', 'o', col_x(2), curr_y); curr_y = curr_y - y_step;
draw_p('$\epsilon$', col_long, '-', 'x', col_x(2), curr_y); curr_y = curr_y - y_step;
draw_p('$k_3$', col_long, '-.', 'none', col_x(2), curr_y);

% --- COLUMN 3: Growth \& Burden ---
curr_y = y_top;
draw_h('Growth \& Burden', 'k', col_x(3), curr_y);
curr_y = curr_y - y_step;
draw_p('$k_1$', col_k1, '-', 'o', col_x(3), curr_y); curr_y = curr_y - y_step;
draw_p('$B_0$', col_B0, '-', 'o', col_x(3), curr_y);

ylim(lg_ax, [0, 1]); xlim(lg_ax, [0, 1]);
end

%% =====================================================================
%% Metric: AUC = int_0^T B(t) dt
%% =====================================================================
function auc = blast_auc(P, IC, Tend, odeOpts)
[t, y] = simulate_beam_ode(P, IC, Tend, odeOpts);
B = y(:,1);
auc = trapz(t, B);
end

%% =====================================================================
%% Deterministic ODE simulation
%% State: [B; E1..EN; A; M]
%% =====================================================================
function [t, y] = simulate_beam_ode(P, IC, Tend, odeOpts)
N = P.N;
y0 = [IC.B0; IC.E0(:); IC.A0; IC.M0];
[t, y] = ode15s(@(tt,yy) beam_ode_full(tt,yy,P), [0 Tend], y0, odeOpts);
y(y < 0) = 0;
end

function dydt = beam_ode_full(~, y, P)
N = P.N;

B = y(1);
E = y(2:1+N);
A = y(2+N);
M = y(3+N);

Etot = sum(E);

% tumour dynamics
dB = P.k1*B*(1 - B/P.K) - P.k2*B*Etot;

% switching fraction
frac = 0;
if (P.Bhalf + B) > 0
    frac = B/(P.Bhalf + B);
end

% effector chain dynamics (robust for N=1,2,...)
dE = zeros(N,1);

if N == 1
    % Single compartment: production into E1 and death out of E1
    dE(1) = 2*P.k4*A*frac - P.delta*E(1);

elseif N == 2
    % Two compartments:
    dE(1) = 2*P.k4*A*frac - P.gamma*E(1);
    dE(2) = 2*P.gamma*E(1) - P.delta*E(2);

else
    % N >= 3 (your original model)
    dE(1) = 2*P.k4*A*frac - P.gamma*E(1);
    for i = 2:N-1
        dE(i) = P.gamma*(2*E(i-1) - E(i));
    end
    dE(N) = 2*P.gamma*E(N-1) - P.delta*E(N);
end

% activated + memory
dA = P.k3*M*B - P.k4*A*(1-frac) - P.k4*A*frac;
dM = -P.k3*M*B + 2*P.k4*A*(1-frac) - P.eps*M;

dydt = [dB; dE; dA; dM];
end

%% =====================================================================
%% Defaults / ICs
%% =====================================================================
function P = default_params()
P.k1    = 0.2;
P.k2    = 2.9e-10;
P.k3    = 1e-9;
P.k4    = 0.1;
P.gamma = 0.3;
P.delta = 0.2;
P.eps   = 0.01;
P.K     = 1e12;
P.N     = 6;
P.Bhalf = 1e9;
end

function IC = default_initial_conditions(P)
CAR0 = 4.1e8;
IC.B0 = 2e11;
IC.A0 = 0;
IC.M0 = 0.468 * CAR0;
IC.E0 = (0.532 * CAR0 / P.N) * ones(P.N,1);
end

function s = prettyDuration(sec)
if sec < 60
    s = sprintf('%.1fs', sec);
elseif sec < 3600
    s = sprintf('%dm %02ds', floor(sec/60), round(mod(sec,60)));
else
    h = floor(sec/3600);
    m = floor(mod(sec,3600)/60);
    ss = round(mod(sec,60));
    s = sprintf('%dh %02dm %02ds', h, m, ss);
end
end