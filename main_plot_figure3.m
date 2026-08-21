function fig3_beam_reproduce()
%FIG3_BEAM_REPRODUCE Reproduce Figure 3 of the BEAM paper (a–d).
%
% Implements BEAM ODE model equations (2.1)–(2.3) and the hybrid method:
% - Deterministic ode45 when B >= Lambda
% - Hybrid when B < Lambda: B simulated as discrete birth/death SSA with
%   propensities depending on continuous states; continuous states integrated
%   deterministically between events with fixed B.

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
k2_grid = linspace(3, 4.5, 151) * 1e-10;

nReps   = 10^3;
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

% (d) v1. This is the version recomputing the data
ax4 = subplot(2,2,4);
plot(ax4, k2_grid/1e-10, Pelim, 'k-', 'LineWidth', 2);
xlabel(ax4, 'CAR T induced tumour death rate ($k_2$)', 'Interpreter','latex');
ylabel(ax4, 'Elimination probability');
xlim(ax4, [3 4.5]);
ylim(ax4, [0 1]);
grid(ax4,'on');
set(ax4, 'Box','on');
title(ax4, 'Elimination probability','FontWeight','normal', 'Color',[0 0 0]);
text(ax4, 0.95, -0.13, '$\times 10^{-10}$', 'Units','normalized','FontSize',14, 'Interpreter','latex');

% (d) v2.
% This is the version with already obtained data
% disp('Careful, using loaded data instead of simulations')
% ax4 = subplot(2,2,4);
% data = load(fullfile(fileparts(mfilename('fullpath')), 'Figure3_data.mat'), 'Pelim', 'k2_grid');
% if isfield(data, 'Pelim')
%     Pelim = data.Pelim;
% end
% k2_grid = linspace(3,4.5, 151);
% plot(ax4, k2_grid, Pelim', 'k-', 'LineWidth', 2);
% xlabel(ax4, 'CAR T induced tumour death rate ($k_2$)', 'Interpreter','latex');
% ylabel(ax4, 'Elimination probability');
% xlim(ax4, [3 4.5]);
% ylim(ax4, [0 1]);
% grid(ax4,'on');
% set(ax4, 'Box','on');
% title(ax4, 'Elimination probability','FontWeight','normal', 'Color',[0 0 0]);
% text(ax4, 0.95, -0.13, '$\times 10^{-10}$', 'Units','normalized','FontSize',14, 'Interpreter','latex');

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
P.k2     = 1e-10;     % overwritten
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
P.A0     = 0;                   
P.E0each = (.532 * P.CAR0)/P.N; % equally split

% Initial blast burden:
P.B0     = 2e11;

end

% =====================================================================
% Hybrid simulation
% State ordering (full ODE): [B; E1..EN; A; M]
% =====================================================================
function out = simulate_beam_hybrid(P, Tend, dt, Lambda)

N  = P.N;

% Time grid for logging (uniform)
tgrid = (0:dt:Tend).';
ng = numel(tgrid);

% Allocate logs
Blog   = zeros(ng,1);
Elog   = zeros(ng,1);
Mlog   = zeros(ng,1);
Alog   = zeros(ng,1);
Tlog   = zeros(ng,1);

% Initial state
B = P.B0;
E = ones(N,1)*P.E0each;
A = P.A0;
M = P.M0;

% Log at t=0
k = 1;
Blog(k) = B;
Elog(k) = sum(E);
Mlog(k) = M;
Alog(k) = A;
Tlog(k) = sum(E) + M + A;

% Run
t = 0;
for k = 2:ng
    tNext = tgrid(k);
    h = tNext - t;

    if B >= Lambda
        % Deterministic step (integrate full system including B)
        y0 = [B; E; A; M];
        opts = odeset('RelTol',1e-7,'AbsTol',1e-10);
        [~,yy] = ode45(@(tt,yy) beam_ode_full(tt,yy,P), [t tNext], y0, opts);
        yend = yy(end,:).';
        B = max(yend(1),0);
        E = max(yend(2:1+N),0);
        A = max(yend(2+N),0);
        M = max(yend(3+N),0);

    else
        % Hybrid step over [t,tNext]:
        % - B treated as discrete birth/death SSA with propensities using current continuous state
        % - Continuous (E,A,M) integrated deterministically between discrete events, holding B fixed
        [B,E,A,M] = hybrid_step_B_discrete(P, B, E, A, M, h);
    end

    t = tNext;

    Blog(k) = B;
    Elog(k) = sum(E);
    Mlog(k) = M;
    Alog(k) = A;
    Tlog(k) = sum(E) + M + A;

    if B <= 0
        % Once eliminated, keep B=0; continue integrating immune contraction deterministically
        B = 0;
    end
end

out.t = tgrid;
out.B = Blog;
out.A = Alog;
out.E = Elog;
out.M = Mlog;      % (M + A) in plot
out.T = Tlog;      % total CAR T = (M + A + sum E_i)
out.eliminated = any(Blog <= 0.5); % elimination = reaches 0 during [0,Tend]
end

% =====================================================================
% One hybrid step when B is discrete
% =====================================================================
function [B,E,A,M] = hybrid_step_B_discrete(P, B, E, A, M, h)

N = P.N;
tloc = 0;

opts = odeset('RelTol',1e-7,'AbsTol',1e-10);

while tloc < h
    % Propensities for discrete B events (using current continuous state)
    Etot = sum(E);
    a_birth = max(P.k1 * B * (1 - B/P.K), 0);
    a_death = max(P.k2 * B * Etot, 0);
    a0 = a_birth + a_death;

    if a0 <= 0
        % No discrete events possible; just integrate continuous over remaining time
        dtc = h - tloc;
        y0 = [E; A; M];
        [~,yy] = ode45(@(tt,yy) beam_ode_continuous(tt,yy,P,B), [0 dtc], y0, opts);
        yend = yy(end,:).';
        E = max(yend(1:N),0);
        A = max(yend(N+1),0);
        M = max(yend(N+2),0);
        tloc = h;
        break;
    end

    % Time to next event
    tau = -log(rand)/a0;

    if tloc + tau > h
        % No event within remaining time; integrate continuous to end
        dtc = h - tloc;
        y0 = [E; A; M];
        [~,yy] = ode45(@(tt,yy) beam_ode_continuous(tt,yy,P,B), [0 dtc], y0, opts);
        yend = yy(end,:).';
        E = max(yend(1:N),0);
        A = max(yend(N+1),0);
        M = max(yend(N+2),0);
        tloc = h;
    else
        % Integrate continuous up to event time
        y0 = [E; A; M];
        [~,yy] = ode45(@(tt,yy) beam_ode_continuous(tt,yy,P,B), [0 tau], y0, opts);
        yend = yy(end,:).';
        E = max(yend(1:N),0);
        A = max(yend(N+1),0);
        M = max(yend(N+2),0);

        % Fire event
        u = rand*a0;
        if u < a_birth
            B = B + 1;
        else
            B = max(B - 1, 0);
        end

        tloc = tloc + tau;
    end
end

end

% =====================================================================
% ODE: full BEAM system (equations (2.1)–(2.3))
% y = [B; E1..EN; A; M]
% =====================================================================
function dydt = beam_ode_full(~, y, P)
N = P.N;

B = y(1);
E = y(2:1+N);
A = y(2+N);
M = y(3+N);

Etot = sum(E);

% (2.1) blasts
dB = P.k1*B*(1 - B/P.K) - P.k2*B*Etot;

% (2.2) effectors
frac = B/(P.Bhalf + B);
dE = zeros(N,1);
dE(1) = 2*P.k4*A*frac - P.gamma*E(1);
for i = 2:N-1
    dE(i) = P.gamma*(2*E(i-1) - E(i));
end
dE(N) = 2*P.gamma*E(N-1) - P.delta*E(N);

% (2.3) active and memory
dA = P.k3*M*B - P.k4*A*(1-frac) - P.k4*A*frac;
dM = -P.k3*M*B + 2*P.k4*A*(1-frac) - P.eps*M;

dydt = [dB; dE; dA; dM];
end

% =====================================================================
% ODE: continuous subsystem when B is held discrete/fixed
% yy = [E1..EN; A; M]
% =====================================================================
function dydt = beam_ode_continuous(~, yy, P, B)
N = P.N;
E = yy(1:N);
A = yy(N+1);
M = yy(N+2);

frac = B/(P.Bhalf + B);

dE = zeros(N,1);
dE(1) = 2*P.k4*A*frac - P.gamma*E(1);
for i = 2:N-1
    dE(i) = P.gamma*(2*E(i-1) - E(i));
end
dE(N) = 2*P.gamma*E(N-1) - P.delta*E(N);

dA = P.k3*M*B - P.k4*A*(1-frac) - P.k4*A*frac;
dM = -P.k3*M*B + 2*P.k4*A*(1-frac) - P.eps*M;

dydt = [dE; dA; dM];
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
