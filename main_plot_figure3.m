function main_plot_figure3()
% Reproduce Figure 3 of the BEAM paper (a–d).
%
% Implements the BEAM reaction network and the reaction-based hybrid method
% described in Supplementary Information, Section A. Each reaction is
% classified independently as deterministic or stochastic; stochastic event
% times are generated from the integrated time-dependent hazard along the
% provisional deterministic trajectory.

% clearvars; close all; clc;
rng(9,'twister');

% =========================
% Parameters (see function below)
% =========================
P = default_params_table1();

% =========================
% Figure 3 settings
% =========================
Tend    = 300;     % days (10 months)
dt      = 0.1;     % hybrid ODE increment 
Lambda  = 1e3;     % switch threshold 
MRD     = 1e6;     % MRD detection threshold
ymin    = 1e3;
ymax    = 1e12;

% k2 values for panels (a,b,c)
k2_vals = [1.5, 2.5, 4] * 1e-10;

% =========================
% Simulate panels (a,b,c)
% =========================
Sims = cell(1,3);
for j = 1:3
    Pj = P; Pj.k2 = k2_vals(j);
    Sims{j} = simulate_beam_hybrid(Pj, Tend, dt, Lambda);
end

% =========================
% Panel (d): elimination probability vs k2
% =========================
k2_grid = linspace(3, 4.5, 2) * 1e-10;

nReps   = 1;
Pelim   = zeros(size(k2_grid));

% Parallel pool if available
p = gcp('nocreate');
if isempty(p)
    try, parpool("threads"); catch, parpool('local'); end
end

fprintf('Running panel (d): %d k2 values × %d realisations...\n', numel(k2_grid), nReps);
q = parallel.pool.DataQueue;
totalJobs = numel(k2_grid)*nReps;
progress = makeProgressPrinter(totalJobs);
afterEach(q, @(~) progress());

for i = 1:numel(k2_grid)
    Ptmp = P; Ptmp.k2 = k2_grid(i);

    elim_here = false(nReps,1);

    parfor r = 1:nReps
        rng(8 + 1000*i + 100000*r, 'twister');
        out = simulate_beam_hybrid(Ptmp, Tend, dt, Lambda);
        elim_here(r) = out.eliminated;
        send(q,1);
    end

    Pelim(i) = mean(elim_here);
end
fprintf('\nDone.\n');

% =========================
% Plot Figure 3 layout (2×2)
% =========================
fig = figure('Color','w','Position',[80 80 1100 760]);

% (a)
ax1 = subplot(2,2,1);
plot_timeseries_panel(ax1, Sims{1}, MRD, ymin, ymax);
title('$k_2 = 1.5 \times 10^{-10}$','FontWeight','normal','Interpreter','latex');
set(ax1,'YScale','log');
yticks(ax1, [1e3 1e6 1e9 1e12]);
set(ax1, 'TickLabelInterpreter', 'tex');
yticklabels(ax1, {'10^{3}','10^{6}','10^{9}','10^{12}'}); 
ylim(ax1,[1e3 1e12]);
lines = findall(ax1, 'Type', 'Line');
set(lines, 'LineWidth', 2);  

% (b)
ax2 = subplot(2,2,2);
plot_timeseries_panel(ax2, Sims{2}, MRD, ymin, ymax);
title('$k_2 = 2.5 \times 10^{-10}$','FontWeight','normal','Interpreter','latex');
set(ax2,'YScale','log');
yticks(ax2, [1e3 1e6 1e9 1e12]);
set(ax2, 'TickLabelInterpreter', 'tex');
yticklabels(ax2, {'10^{3}','10^{6}','10^{9}','10^{12}'});
ylim(ax2,[1e3 1e12]);
lines = findall(ax2, 'Type', 'Line');
set(lines, 'LineWidth', 2); 

% (c)
ax3 = subplot(2,2,3);
plot_timeseries_panel(ax3, Sims{3}, MRD, ymin, ymax);
title('$k_2 = 4 \times 10^{-10}$','FontWeight','normal','Interpreter','latex');
set(ax3,'YScale','log');
yticks(ax3, [1e3 1e6 1e9 1e12]);
set(ax3, 'TickLabelInterpreter', 'tex');
yticklabels(ax3, {'10^{3}','10^{6}','10^{9}','10^{12}'});
ylim(ax3,[1e3 1e12]);
lines = findall(ax3, 'Type', 'Line');
set(lines, 'LineWidth', 2);

% (d) This is the version recomputing the data
% plot(ax4, k2_grid/1e-10, Pelim, 'k-', 'LineWidth', 2);
% xlabel(ax4, 'CAR T induced tumour death rate ($k_2$)', 'Interpreter','latex');
% ylabel(ax4, 'Elimination probability');
% xlim(ax4, [3 4.5]);
% ylim(ax4, [0 1]);
% grid(ax4,'on');
% set(ax4, 'Box','on');
% title(ax4, 'Elimination probability','FontWeight','normal', 'Color',[0 0 0]);
% text(ax4, 0.95, -0.13, '$\times 10^{-10}$', 'Units','normalized','FontSize',14, 'Interpreter','latex');

% (d) With pre-computed data
disp('Using loaded data instead of simulations')
ax4 = subplot(2,2,4);
data = load(fullfile(fileparts(mfilename('fullpath')), 'Figure3_data.mat'), 'Pelim', 'k2_grid');
if isfield(data, 'Pelim')
    Pelim = data.Pelim;
end
k2_grid = linspace(3,4.5, 151);
plot(ax4, k2_grid, Pelim', 'k-', 'LineWidth', 2);
xlabel(ax4, 'CAR T induced tumour death rate ($k_2$)', 'Interpreter','latex');
ylabel(ax4, 'Elimination probability');
xlim(ax4, [3 4.5]);
ylim(ax4, [0 1]);
grid(ax4,'on');
set(ax4, 'Box','on');
title(ax4, 'Elimination probability','FontWeight','normal', 'Color',[0 0 0]);
text(ax4, 0.95, -0.13, '$\times 10^{-10}$', 'Units','normalized','FontSize',14, 'Interpreter','latex');

% Global styling
set([ax1 ax2 ax3 ax4], 'FontSize', 11, 'LineWidth', 1);
text(ax1, -0.16, 1.05, '(a)', 'Units','normalized','FontSize',14);
text(ax2, -0.16, 1.05, '(b)', 'Units','normalized','FontSize',14);
text(ax3, -0.16, 1.05, '(c)', 'Units','normalized','FontSize',14);
text(ax4, -0.16, 1.05, '(d)', 'Units','normalized','FontSize',14);
% Save
% outbase = sprintf('Figure3_BEAM_%s', datestr(now,'yyyymmdd_HHMMSS'));
% print(fig, [outbase '.png'], '-dpng', '-r300');
% print(fig, [outbase '.pdf'], '-dpdf', '-vector');
% save([outbase '.mat'], 'P','Sims','k2_grid','Pelim','MRD','Lambda','dt','Tend');
% savefig(fig, [outbase '.fig']);
% 
% fprintf('Saved figure/data as %s.{png,pdf,mat}\n', outbase);

end

% =====================================================================
% Defaults: Table 1 baseline parameters + initialisation logic
% =====================================================================
function P = default_params_table1()
% Baseline parameters (Table 1) and text (initialisation fractions).
P.k1     = 0.2;
P.k2     = 4e-10;     % Table 1 baseline; overwritten where required
P.k3     = 1e-9;
P.k4     = 0.1;
P.gamma  = 0.3;
P.delta  = 0.2;
P.eps    = 0.01;
P.K      = 1e12;
P.N      = 6;

P.Bhalf  = 1e9;

% Initial conditions:
% CAR T total dose = 4.1e8 at t=0
% Memory M gets naive+central memory = 31.2% + 15.6% = 46.8%
% Effector compartments E1..E6 share effector memory + terminal effector = 19.2% + 34.0% = 53.2% equally.
P.CAR0   = 4.1e8;
P.M0     = .468 * P.CAR0;
P.A0     = 0e4;                   
P.E0each = (.532 * P.CAR0)/P.N; % equally split

% Initial blast burden:
P.B0     = 2e11;

end

% =====================================================================
% Reaction-based hybrid simulation (Supplementary Information, Section A)
% State ordering: x = [B; E1; ...; EN; A; M]
% =====================================================================
function out = simulate_beam_hybrid(P, Tend, dt, Lambda)

N = P.N;

% Uniform output / macro-step grid. In this reproduction code we use
% dt = 0.1 day both as the macro-step size Delta t and as the fast-reaction
% classification threshold delta t appearing in the Supplement.
tgrid = (0:dt:Tend).';
ng = numel(tgrid);

% State indices
idxB = 1;
idxE = 2:(N+1);
idxA = N + 2;
idxM = N + 3;
nState = N + 3;

% Stoichiometric matrix for the N+7 reactions in Supplementary Table 1
S = beam_stoichiometry(P);

% Initial state
x = zeros(nState,1);
x(idxB) = P.B0;
x(idxE) = P.E0each;
x(idxA) = P.A0;
x(idxM) = P.M0;

% Allocate output logs
Blog = zeros(ng,1);
Elog = zeros(ng,1);
Mlog = zeros(ng,1);
Alog = zeros(ng,1);
Tlog = zeros(ng,1);

% Log initial state
Blog(1) = x(idxB);
Elog(1) = sum(x(idxE));
Mlog(1) = x(idxM);
Alog(1) = x(idxA);
Tlog(1) = sum(x(idxE)) + x(idxA) + x(idxM);

% ODE tolerances used for the deterministic reaction subset.
% NonNegative prevents small negative values generated by ODE tolerances.
odeOpts = odeset('RelTol',1e-7,'AbsTol',1e-10, ...
                  'NonNegative',1:nState);

for k = 2:ng
    t0 = tgrid(k-1);
    t1 = tgrid(k);

    % Advance one macro-step. Reactions are classified at t0 and then
    % reclassified after every stochastic event occurring inside [t0,t1].
    x = hybrid_reaction_interval(P, x, t0, t1, dt, Lambda, S, odeOpts);

    % Numerical safeguard for tiny negative round-off values.
    x = max(x,0);

    Blog(k) = x(idxB);
    Elog(k) = sum(x(idxE));
    Mlog(k) = x(idxM);
    Alog(k) = x(idxA);
    Tlog(k) = sum(x(idxE)) + x(idxA) + x(idxM);
end

out.t = tgrid;
out.B = Blog;
out.A = Alog;
out.E = Elog;
out.M = Mlog;
out.T = Tlog;
out.eliminated = any(Blog == 0);
end

% =====================================================================
% Advance one macro-step of the reaction-based hybrid algorithm
% =====================================================================
function x = hybrid_reaction_interval(P, x, tStart, tEnd, dtClass, Lambda, S, odeOpts)

nState = numel(x);
t = tStart;

while t < tEnd
    % ---------------------------------------------------------------
    % 1. Classify every reaction from the state at the current restart
    %    time. A reaction is deterministic iff
    %       (i) 1/alpha_j < dtClass, and
    %      (ii) every reactant population exceeds Lambda.
    %    All other reactions are stochastic.
    % ---------------------------------------------------------------
    alpha0 = beam_propensities(x,P);
    detMask = classify_reactions(x, alpha0, P, dtClass, Lambda);
    stochMask = ~detMask;
     if x(1) > 0 && x(1) <= Lambda
        x(1) = stochastic_round(x(1));
    end

    % If no stochastic reaction is present, integrate deterministic
    % reactions directly to the end of the macro-step.
    if ~any(stochMask)
        [~,xx] = ode45(@(tt,xx) deterministic_rhs(tt,xx,P,S,detMask), ...
                       [t tEnd], x, odeOpts);
        x = xx(end,:).';
        t = tEnd;
        continue;
    end

    % ---------------------------------------------------------------
    % 2. Draw the unit-exponential threshold and integrate the
    %    provisional deterministic trajectory together with the
    %    accumulated stochastic hazard
    %
    %        H'(s) = sum_{j in S(t)} alpha_j(x(s)).
    %
    %    The next event occurs when H = -log(1-u).
    % ---------------------------------------------------------------
    u = rand;
    hazardTarget = -log(max(1-u, realmin));

    z0 = [x; 0];
    eventOpts = odeset(odeOpts, ...
        'Events', @(tt,zz) hazard_event(tt,zz,hazardTarget));

    [~,zz,te,ze] = ode45( ...
        @(tt,zz) augmented_hybrid_rhs(tt,zz,P,S,detMask,stochMask,nState), ...
        [t tEnd], z0, eventOpts);

    if isempty(te)
        % No stochastic event before the end of this macro-step.
        x = zz(end,1:nState).';
        t = tEnd;
        continue;
    end

    % ---------------------------------------------------------------
    % 3. A stochastic event occurred. Evaluate the stochastic
    %    propensities at the event time and choose reaction i with
    %
    %        P(i | event at t+tau) = alpha_i / sum_j alpha_j.
    % ---------------------------------------------------------------
    xEvent = ze(end,1:nState).';
    alphaEvent = beam_propensities(xEvent,P);
    rates = alphaEvent;
    rates(~stochMask) = 0;
    totalRate = sum(rates);

    if totalRate <= 0
        % This should only arise from numerical event-location round-off.
        % Advance to the located event time and reclassify without a jump.
        x = xEvent;
        t = te(end);
        continue;
    end

    r = rand * totalRate;
    j = find(cumsum(rates) >= r, 1, 'first');

    % Apply the stoichiometric jump exactly and restart/reclassify.
    x = xEvent + S(:,j);

    % A deterministic component may be fractional when it is reclassified
    % as stochastic. If a final -1 jump crosses zero by a sub-cell amount,
    % interpret this as extinction and prevent an unphysical negative count.
    x = max(x,0);

    t = te(end);
end
end

% =====================================================================
% Deterministic contribution from the currently deterministic reactions
% =====================================================================
function dx = deterministic_rhs(~, x, P, S, detMask)
alpha = beam_propensities(x,P);
if any(detMask)
    dx = S(:,detMask) * alpha(detMask);
else
    dx = zeros(size(x));
end
end

% =====================================================================
% Provisional deterministic state + integrated stochastic hazard
% =====================================================================
function dz = augmented_hybrid_rhs(~, z, P, S, detMask, stochMask, nState)
x = z(1:nState);
alpha = beam_propensities(x,P);

if any(detMask)
    dx = S(:,detMask) * alpha(detMask);
else
    dx = zeros(nState,1);
end

dH = sum(alpha(stochMask));
dz = [dx; dH];
end

% =====================================================================
% Event: accumulated stochastic hazard reaches the exponential threshold
% =====================================================================
function [value,isterminal,direction] = hazard_event(~, z, hazardTarget)
value = z(end) - hazardTarget;
isterminal = 1;
direction = 1;
end

% =====================================================================
% Reaction propensities from Supplementary Table 1
% Reactions:
% R1      B -> B+1                          k1 B
% R2      B -> B-1                          k1 B^2 / K
% R3      B -> B-1                          k2 B E
% R4      M -> M-1, A -> A+1                k3 M B
% R5      A -> A-1, M -> M+2                k4 A Bhalf/(Bhalf+B)
% R6      A -> A-1, E1 -> E1+2              k4 A B/(Bhalf+B)
% R6+i    Ei -> Ei-1, E(i+1) -> E(i+1)+2    gamma Ei
% RN+6    EN -> EN-1                        delta EN
% RN+7    M -> M-1                          eps M
% =====================================================================
function alpha = beam_propensities(x,P)

N = P.N;
x = max(x,0);  % protects propensity evaluation from ODE round-off

B = x(1);
E = x(2:N+1);
A = x(N+2);
M = x(N+3);
Etot = sum(E);

fracEff = B/(P.Bhalf + B);
fracMem = P.Bhalf/(P.Bhalf + B);

alpha = zeros(N+7,1);
alpha(1) = P.k1 * B;
alpha(2) = P.k1 * B^2 / P.K;
alpha(3) = P.k2 * B * Etot;
alpha(4) = P.k3 * M * B;
alpha(5) = P.k4 * A * fracMem;
alpha(6) = P.k4 * A * fracEff;

for i = 1:N-1
    alpha(6+i) = P.gamma * E(i);
end

alpha(N+6) = P.delta * E(N);
alpha(N+7) = P.eps * M;

% Propensities must be finite and non-negative.
alpha(~isfinite(alpha) | alpha < 0) = 0;
end

% =====================================================================
% Stoichiometric update matrix corresponding to Supplementary Table 1
% =====================================================================
function S = beam_stoichiometry(P)

N = P.N;
nState = N + 3;
nRxn = N + 7;
S = zeros(nState,nRxn);

idxB = 1;
idxA = N + 2;
idxM = N + 3;

% R1: blast birth
S(idxB,1) = 1;

% R2: density-dependent blast death
S(idxB,2) = -1;

% R3: CAR T-cell-induced blast death
S(idxB,3) = -1;

% R4: memory activation, M -> A
S(idxM,4) = -1;
S(idxA,4) = 1;

% R5: activated cell divides into two memory cells
S(idxA,5) = -1;
S(idxM,5) = 2;

% R6: activated cell divides into two E1 cells
S(idxA,6) = -1;
S(2,6) = 2;

% R6+i: Ei divides into two E(i+1), i = 1,...,N-1
for i = 1:N-1
    j = 6 + i;
    idxEi = 1 + i;
    idxEip1 = 2 + i;
    S(idxEi,j) = -1;
    S(idxEip1,j) = 2;
end

% RN+6: terminal effector death
S(N+1,N+6) = -1;

% RN+7: memory death
S(idxM,N+7) = -1;
end


% =====================================================================
% Unbiased stochastic rounding to a non-negative integer
% =====================================================================
function n = stochastic_round(x)

if x <= 0
    n = 0;
    return;
end

n0 = floor(x);
fraction = x - n0;

if rand < fraction
    n = n0 + 1;
else
    n = n0;
end

end


% =====================================================================
% Reaction-by-reaction deterministic/stochastic partition
% =====================================================================
function detMask = classify_reactions(x, alpha, P, dtClass, Lambda)

N = P.N;
B = max(x(1),0);
E = max(x(2:N+1),0);
A = max(x(N+2),0);
M = max(x(N+3),0);
Etot = sum(E);

% Criterion (i): expected next reaction time 1/alpha_j < dtClass.
fast = false(N+7,1);
pos = alpha > 0;
fast(pos) = (1 ./ alpha(pos)) < dtClass;

% Criterion (ii): every reactant population exceeds Lambda.
reactantsLarge = false(N+7,1);
reactantsLarge(1) = B > Lambda;                     % R1: B
reactantsLarge(2) = B > Lambda;                     % R2: B
reactantsLarge(3) = B > Lambda && Etot > Lambda;    % R3: B and total E
reactantsLarge(4) = B > Lambda && M > Lambda;       % R4: B and M
reactantsLarge(5) = A > Lambda;                     % R5: A
reactantsLarge(6) = A > Lambda;                     % R6: A

for i = 1:N-1
    reactantsLarge(6+i) = E(i) > Lambda;             % R6+i: Ei
end

reactantsLarge(N+6) = E(N) > Lambda;                % terminal E_N death
reactantsLarge(N+7) = M > Lambda;                   % memory death

detMask = fast & reactantsLarge;
end


% =====================================================================
% Plot helper: time-series panel matching Figure 3 styling
% =====================================================================
function plot_timeseries_panel(ax, sim, MRD, ymin, ymax)

axes(ax);
hold on;

% Background MRD bands
x0 = sim.t(1); x1 = sim.t(end);
% Dark band: MRD negative (<= MRD)
p1 = patch([x0 x1 x1 x0], [ymin ymin MRD MRD], [0.85 0.85 0.85], ...
    'EdgeColor','none');
% Light band: MRD positive (MRD to 1e10) — matches the look of the figure
p2 = patch([x0 x1 x1 x0], [MRD MRD 1e10 1e10], [0.92 0.92 0.92], ...
    'EdgeColor','none');
% push shading to the very bottom
uistack(p1,'bottom'); 
uistack(p2,'bottom');
% do not include in the legend
p1.HandleVisibility = 'off';
p2.HandleVisibility = 'off';
% make sure grid is above (MATLAB default is usually top, but force it)
ax.Layer = 'top';

% Curves
plot(sim.t, sim.B, 'r-', 'LineWidth', 1.3);
plot(sim.t, sim.E, 'Color',[0 0.8 0], 'LineWidth', 1.3); % green
plot(sim.t, sim.M + sim.A, 'b-', 'LineWidth', 1.3);
plot(sim.t, sim.T, 'k--', 'LineWidth', 1.2);

set(ax, 'YScale','log');
xlim(ax, [0 300]);
ylim(ax, [ymin ymax]);

xlabel(ax, 'Time [days]');
ylabel(ax, 'Cell Number');

grid(ax,'on');
set(ax,'Box','on');

text(200, 2.5e6, 'MRD positive', 'FontSize', 9, 'Color',[0.3 0.3 0.3]);
text(200, 6e5,  'MRD negative', 'FontSize', 9, 'Color',[0.3 0.3 0.3]);

legend({'Blast','Effector','Memory','Total'}, 'Location','northeast', 'Box','on', 'NumColumns', 2);
hold off;
end

% =====================================================================
% Progress printer
% =====================================================================
function update = makeProgressPrinter(total)
startTime = tic; counter = 0;
    function cb()
        counter = counter + 1;
        pct = 100*counter/total;
        elapsed = toc(startTime);
        rate = counter / max(elapsed, eps);
        remaining = (total - counter) / max(rate, eps);
        fprintf('\rProgress: %d/%d (%.1f%%) | Elapsed: %s | ETA: %s', ...
            counter, total, pct, prettyDuration(elapsed), prettyDuration(remaining));
        if counter == total, fprintf('\n'); end
    end
update = @cb;
end

function s = prettyDuration(sec)
if sec < 60
    s = sprintf('%.1fs', sec);
elseif sec < 3600
    s = sprintf('%dm %02ds', floor(sec/60), round(mod(sec,60)));
else
    h = floor(sec/3600); m = floor(mod(sec,3600)/60); ssec = round(mod(sec,60));
    s = sprintf('%dh %02dm %02ds', h, m, ssec);
end
end
