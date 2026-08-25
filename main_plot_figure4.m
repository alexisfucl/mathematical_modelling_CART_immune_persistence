function fig4_beam()
% Reproduces Figure 4 of the BEAM paper
% Implements the mathematically exact Integral Method for the PDMP hybrid
% model. 

rng(9,'twister');

% =========================
% Parameters (Table 1 baseline)
% =========================
P = default_params_table1();

Tend    = 1000;     
dt      = 0.1;      
Lambda  = 1e3;     
MRD     = 1e6;     
ymin    = 1e0;
ymax    = 1e12;

k1_vals   = [0.06, 0.1, 0.16];

% =========================
% Simulate panels (a,b,c)
% =========================
Sims = cell(1,3);
for j = 1:3
    disp('start')
    Pj = P; Pj.k1 = k1_vals(j);
    Sims{j} = simulate_beam_exact_pdmp(Pj, Tend, dt, Lambda);
    disp('end')
end

% =========================
% Plot Figure 4 layout (1×3)
% =========================
fig = figure('Color','w','Position',[1500 1500 200 200]);
t = tiledlayout(1,3); 

t.TileSpacing = 'compact'; 
t.Padding = 'tight';
position=1;
ax1 = subplot(1,3,1); plot_timeseries_panel(ax1, Sims{1}, MRD, ymin, ymax, position);
title('$k_1 = 0.06$','FontWeight','normal','Interpreter','latex');
set(ax1,'YScale','log'); ylim(ax1,[1e0 1e12]);
yticks(ax1, [1e0 1e3 1e6 1e9 1e12]); 
set(ax1, 'TickLabelInterpreter', 'tex');
yticklabels(ax1, {'10^0','10^3','10^6','10^9','10^{12}'});  
position=2;
ax2 = subplot(1,3,2); plot_timeseries_panel(ax2, Sims{2}, MRD, ymin, ymax, position);
title('$k_1 = 0.1$','FontWeight','normal','Interpreter','latex');
set(ax2,'YScale','log'); ylim(ax2,[1e0 1e12]);
yticks(ax2, [1e0 1e3 1e6 1e9 1e12]); 
set(ax2, 'TickLabelInterpreter', 'tex');
yticklabels(ax2, {'10^0','10^3','10^6','10^9','10^{12}'});  
position=1;
ax3 = subplot(1,3,3); plot_timeseries_panel(ax3, Sims{3}, MRD, ymin, ymax, position);
title('$k_1 = 0.16$','FontWeight','normal','Interpreter','latex');
set(ax3,'YScale','log'); ylim(ax3,[1e0 1e12]);
yticks(ax3, [1e0 1e3 1e6 1e9 1e12]); 
set(ax3, 'TickLabelInterpreter', 'tex');
yticklabels(ax3, {'10^0','10^3','10^6','10^9','10^{12}'});  

set([ax1 ax2 ax3], 'FontSize', 11, 'LineWidth', 1);
lines = findall(fig, 'Type', 'Line'); set(lines, 'LineWidth', 2);
text(ax1, -0.16, 1.05, '(a)', 'Units','normalized','FontSize',14);
text(ax2, -0.16, 1.05, '(b)', 'Units','normalized','FontSize',14);
text(ax3, -0.16, 1.05, '(c)', 'Units','normalized','FontSize',14);

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 20, 4]); 
set(findall(gcf,'type','axes'), 'PlotBoxAspectRatioMode', 'auto');
exportgraphics(gcf, 'Figure4.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
end

% =====================================================================
% Defaults
% =====================================================================
function P = default_params_table1()
P.k1     = 0.2;     
P.k2     = 2.9e-10; 
P.k3     = 1e-9;
P.k4     = 0.1;
P.gamma  = 0.3;
P.delta  = 0.2;
P.eps    = 0.014;
P.K      = 1e12;
P.N      = 6;
P.Bhalf  = 1e9;  

P.CAR0   = 4.1e8;
P.M0     = .468 * P.CAR0;
P.A0     = 0;                  
P.E0each = (.532 * P.CAR0)/P.N; 
P.B0     = 2e11;
end

% =====================================================================
% Reaction-based hybrid simulation (Supplementary Information, Section A)
% State ordering: x = [B; E1; ...; EN; A; M]
% =====================================================================
function out = simulate_beam_exact_pdmp(P, Tend, dt, Lambda)

N = P.N;
tgrid = (0:dt:Tend).';
ng = numel(tgrid);

idxB = 1;
idxE = 2:(N+1);
idxA = N + 2;
idxM = N + 3;
nState = N + 3;

% Stoichiometric matrix for the N+7 reactions in Supplementary Table 1.
S = beam_stoichiometry(P);

% Initial state.
x = zeros(nState,1);
x(idxB) = P.B0;
x(idxE) = P.E0each;
x(idxA) = P.A0;
x(idxM) = P.M0;

% Track whether each species was fully stochastic at the previous
% classification point. This is used only for the continuous-to-discrete
% hand-off, so that a newly fully stochastic species is integerised once.
wasFullyStochastic = false(nState,1);

% Output logs.
Blog = zeros(ng,1);
Elog = zeros(ng,1);
Mlog = zeros(ng,1);
Tlog = zeros(ng,1);

Blog(1) = x(idxB);
Elog(1) = sum(x(idxE));
Mlog(1) = x(idxM) + x(idxA);
Tlog(1) = sum(x(idxE)) + x(idxA) + x(idxM);

% Deterministic reactions are integrated with ode45, as described in the
% Supplementary Information.
odeOpts = odeset('RelTol',1e-7,'AbsTol',1e-10, ...
                  'NonNegative',1:nState);

for k = 2:ng
    t0 = tgrid(k-1);
    t1 = tgrid(k);

    [x,wasFullyStochastic] = hybrid_reaction_interval( ...
        P, x, t0, t1, dt, Lambda, S, odeOpts, wasFullyStochastic);

    x = max(x,0);

    Blog(k) = x(idxB);
    Elog(k) = sum(x(idxE));
    Mlog(k) = x(idxM) + x(idxA);
    Tlog(k) = sum(x(idxE)) + x(idxA) + x(idxM);
end

out.t = tgrid;
out.B = Blog;
out.E = Elog;
out.M = Mlog;
out.T = Tlog;
out.eliminated = any(Blog == 0);
end


% =====================================================================
% Advance one macro-step of the reaction-based hybrid algorithm
% =====================================================================
function [x,wasFullyStochastic] = hybrid_reaction_interval( ...
    P, x, tStart, tEnd, dtClass, Lambda, S, odeOpts, wasFullyStochastic)

nState = numel(x);
t = tStart;
timeTol = 100*eps(max(1,abs(tEnd)));

while t < tEnd - timeTol

    % -------------------------------------------------------------
    % 1. Classify every reaction independently.
    % A reaction is deterministic only when BOTH
    %   (i)  1/alpha_j < dtClass, and
    %   (ii) all of its reactant populations exceed Lambda.
    % -------------------------------------------------------------
    alpha = beam_propensities(x,P);
    detMask = classify_reactions(x,alpha,P,dtClass,Lambda);
    stochMask = ~detMask;

    % -------------------------------------------------------------
    % Continuous-to-discrete hand-off.
    % -------------------------------------------------------------
    fullyStochastic = ~any(abs(S(:,detMask)) > 0,2);
    newlyFullyStochastic = fullyStochastic & ~wasFullyStochastic;

    if any(newlyFullyStochastic)
        ids = find(newlyFullyStochastic);
        for q = 1:numel(ids)
            x(ids(q)) = stochastic_round_nonnegative(x(ids(q)));
        end

        % Rounding can very slightly change a propensity, so reclassify.
        alpha = beam_propensities(x,P);
        detMask = classify_reactions(x,alpha,P,dtClass,Lambda);
        stochMask = ~detMask;
        fullyStochastic = ~any(abs(S(:,detMask)) > 0,2);
    end
    wasFullyStochastic = fullyStochastic;

    % -------------------------------------------------------------
    % 2. If every reaction is deterministic, integrate directly to
    %    the end of the current macro-step.
    % -------------------------------------------------------------
    if ~any(stochMask)
        [~,xx] = ode45(@(tt,xx) deterministic_rhs(tt,xx,P,S,detMask), ...
                       [t tEnd], x, odeOpts);
        x = xx(end,:).';
        t = tEnd;
        continue;
    end

    % -------------------------------------------------------------
    % 3. Draw the exponential threshold for the next stochastic event.
    % -------------------------------------------------------------
    u = rand;
    hazardTarget = -log(max(1-u,realmin));

    % Append the accumulated stochastic hazard H to the state.
    z0 = [x; 0];
    eventOpts = odeset(odeOpts,'Events', ...
        @(tt,zz) hazard_event(tt,zz,hazardTarget));

    % The reaction partition is held fixed until either a stochastic event
    % occurs or the current macro-step ends, exactly as described in the
    % Supplementary Information.
    [~,zz,te,ze] = ode45( ...
        @(tt,zz) augmented_hybrid_rhs(tt,zz,P,S,detMask,stochMask,nState), ...
        [t tEnd], z0, eventOpts);

    if isempty(te)
        % The integrated stochastic hazard did not reach the threshold.
        x = zz(end,1:nState).';
        t = tEnd;
        continue;
    end

    % -------------------------------------------------------------
    % 4. A stochastic reaction fires at t = te(end).
    %    Re-evaluate propensities at the event state and choose reaction j
    %    with probability alpha_j/sum(alpha_i) over the stochastic set.
    % -------------------------------------------------------------
    xEvent = ze(end,1:nState).';
    alphaEvent = beam_propensities(xEvent,P);
    rates = feasible_stochastic_rates(xEvent,alphaEvent,S,stochMask);
    totalRate = sum(rates);

    if totalRate <= 0
        % Numerical safeguard: if the stochastic rate has vanished exactly
        % at the located event, restart from that state and reclassify.
        x = xEvent;
        t = te(end);
        continue;
    end

    r = rand * totalRate;
    j = find(cumsum(rates) >= r,1,'first');

    % Apply the exact stoichiometric jump.
    x = xEvent + S(:,j);

    % Guard only against tiny numerical undershoots in mixed continuous /
    % stochastic species. Fully discrete species remain integer-valued.
    x(x < 0 & x > -1e-9) = 0;
    if any(x < -1e-9)
        error('Hybrid update produced a negative population after reaction %d.',j);
    end

    t = te(end);
    % Loop back immediately: all reactions are reclassified after the jump.
end
end


% =====================================================================
% Deterministic contribution from the currently deterministic reactions
% =====================================================================
function dx = deterministic_rhs(~,x,P,S,detMask)
alpha = beam_propensities(x,P);
if any(detMask)
    dx = S(:,detMask) * alpha(detMask);
else
    dx = zeros(size(x));
end
end


% =====================================================================
% Deterministic state plus accumulated stochastic hazard
% =====================================================================
function dz = augmented_hybrid_rhs(~,z,P,S,detMask,stochMask,nState)
x = z(1:nState);
alpha = beam_propensities(x,P);

if any(detMask)
    dx = S(:,detMask) * alpha(detMask);
else
    dx = zeros(nState,1);
end

rates = feasible_stochastic_rates(x,alpha,S,stochMask);
dH = sum(rates);
dz = [dx; dH];
end


% =====================================================================
% Positivity safeguard for stochastic reactions in mixed species
% =====================================================================
function rates = feasible_stochastic_rates(x,alpha,S,stochMask)
rates = alpha;
rates(~stochMask) = 0;

for j = find(stochMask).'
    required = max(-S(:,j),0);
    if any(x + 1e-12 < required)
        rates(j) = 0;
    end
end
end


% =====================================================================
% Event condition: stop when the integrated stochastic hazard reaches R
% =====================================================================
function [value,isterminal,direction] = hazard_event(~,z,hazardTarget)
value = z(end) - hazardTarget;
isterminal = 1;
direction = 1;
end


% =====================================================================
% Reaction propensities from Supplementary Table 1
% =====================================================================
function alpha = beam_propensities(x,P)
N = P.N;
x = max(x,0);

B = x(1);
E = x(2:N+1);
A = x(N+2);
M = x(N+3);
Etot = sum(E);

if P.Bhalf + B > 0
    fracEff = B/(P.Bhalf + B);
    fracMem = P.Bhalf/(P.Bhalf + B);
else
    fracEff = 0;
    fracMem = 1;
end

alpha = zeros(N+7,1);

alpha(1) = P.k1 * B;
alpha(2) = P.k1 * B^2/P.K;
alpha(3) = P.k2 * B * Etot;
alpha(4) = P.k3 * M * B;
alpha(5) = P.k4 * A * fracMem;
alpha(6) = P.k4 * A * fracEff;

for i = 1:N-1
    alpha(6+i) = P.gamma * E(i);
end

alpha(N+6) = P.delta * E(N);
alpha(N+7) = P.eps * M;

alpha(~isfinite(alpha) | alpha < 0) = 0;
end


% =====================================================================
% Stoichiometric matrix corresponding to Supplementary Table 1
% =====================================================================
function S = beam_stoichiometry(P)
N = P.N;
nState = N + 3;
nRxn = N + 7;

S = zeros(nState,nRxn);

idxB = 1;
idxA = N + 2;
idxM = N + 3;

% R1: B -> B + 1
S(idxB,1) = 1;

% R2: B -> B - 1
S(idxB,2) = -1;

% R3: B -> B - 1
S(idxB,3) = -1;

% R4: M -> M - 1, A -> A + 1
S(idxM,4) = -1;
S(idxA,4) = 1;

% R5: A -> A - 1, M -> M + 2
S(idxA,5) = -1;
S(idxM,5) = 2;

% R6: A -> A - 1, E1 -> E1 + 2
S(idxA,6) = -1;
S(2,6) = 2;

% R6+i: Ei -> Ei - 1, E(i+1) -> E(i+1) + 2
for i = 1:N-1
    j = 6 + i;
    idxEi = 1 + i;
    idxEip1 = 2 + i;
    S(idxEi,j) = -1;
    S(idxEip1,j) = 2;
end

% R(N+6): EN -> EN - 1
S(N+1,N+6) = -1;

% R(N+7): M -> M - 1
S(idxM,N+7) = -1;
end


% =====================================================================
% Reaction-by-reaction deterministic/stochastic classification
% =====================================================================
function detMask = classify_reactions(x,alpha,P,dtClass,Lambda)
N = P.N;

B = max(x(1),0);
E = max(x(2:N+1),0);
A = max(x(N+2),0);
M = max(x(N+3),0);
Etot = sum(E);

% Criterion 1: expected next-event time 1/alpha_j < dtClass.
fast = false(N+7,1);
pos = alpha > 0;
fast(pos) = (1./alpha(pos)) < dtClass;

% Criterion 2: all reactant populations exceed Lambda.
reactantsLarge = false(N+7,1);
reactantsLarge(1) = B > Lambda;                    % R1
reactantsLarge(2) = B > Lambda;                    % R2
reactantsLarge(3) = B > Lambda && Etot > Lambda;   % R3
reactantsLarge(4) = B > Lambda && M > Lambda;      % R4
reactantsLarge(5) = A > Lambda;                    % R5
reactantsLarge(6) = A > Lambda;                    % R6

for i = 1:N-1
    reactantsLarge(6+i) = E(i) > Lambda;
end

reactantsLarge(N+6) = E(N) > Lambda;
reactantsLarge(N+7) = M > Lambda;

detMask = fast & reactantsLarge;
end


% =====================================================================
% Unbiased stochastic rounding for a newly fully stochastic species
% =====================================================================
function n = stochastic_round_nonnegative(x)
if x <= 0
    n = 0;
    return;
end

n0 = floor(x);
p = x - n0;
n = n0 + (rand < p);
end

% =====================================================================
% Plot helper
% =====================================================================
function plot_timeseries_panel(ax, sim, MRD, ymin, ymax, position)
axes(ax); %#ok<LAXES>
hold on;
x0 = sim.t(1); x1 = sim.t(end);
p1 = patch([x0 x1 x1 x0], [ymin ymin MRD MRD], [0.85 0.85 0.85], 'EdgeColor','none');
p2 = patch([x0 x1 x1 x0], [MRD MRD 1e10 1e10], [0.92 0.92 0.92], 'EdgeColor','none');
uistack(p1,'bottom'); uistack(p2,'bottom');
p1.HandleVisibility = 'off'; p2.HandleVisibility = 'off';
ax.Layer = 'top';

plot(sim.t, sim.B, 'r-', 'LineWidth', 1.3);
plot(sim.t, sim.E, 'Color',[0 0.8 0], 'LineWidth', 1.3); 
plot(sim.t, sim.M, 'b-', 'LineWidth', 1.3);
plot(sim.t, sim.T, 'k--', 'LineWidth', 1.2);

set(ax, 'YScale','log'); xlim(ax, [0 1000]); ylim(ax, [ymin ymax]);
xlabel(ax, 'Time [days]'); ylabel(ax, 'Cell Number');
grid(ax,'on'); set(ax,'Box','on');
if position==1
    text(750, 2.5e6, 'MRD positive', 'FontSize', 9, 'Color',[0.3 0.3 0.3]);
    text(750, 6e5,  'MRD negative', 'FontSize', 9, 'Color',[0.3 0.3 0.3]);
else
    text(250, 2.5e6, 'MRD positive', 'FontSize', 9, 'Color',[0.3 0.3 0.3]);
    text(250, 6e9,  'MRD negative', 'FontSize', 9, 'Color',[0.3 0.3 0.3]);
end
legend({'Blast','Effector','Memory','Total'}, 'Location','northeast', 'Box','on', 'NumColumns', 2);
hold off;
end

% =====================================================================
% Progress printer
% =====================================================================
function update = makeProgressPrinter(total)
startTime = tic; counter = 0;
    function cb()
        counter = counter + 1; pct = 100*counter/total; elapsed = toc(startTime);
        rate = counter / max(elapsed, eps); remaining = (total - counter) / max(rate, eps);
        fprintf('\rProgress: %d/%d (%.1f%%) | Elapsed: %s | ETA: %s', ...
            counter, total, pct, prettyDuration(elapsed), prettyDuration(remaining));
        if counter == total, fprintf('\n'); end
    end
update = @cb;
end

function s = prettyDuration(sec)
if sec < 60, s = sprintf('%.1fs', sec);
elseif sec < 3600, s = sprintf('%dm %02ds', floor(sec/60), round(mod(sec,60)));
else, h = floor(sec/3600); m = floor(mod(sec,3600)/60); ssec = round(mod(sec,60)); s = sprintf('%dh %02dm %02ds', h, m, ssec);
end
end