clearvars; close all; clc;

% Create a large wide figure for the 3-panel layout
figure('Color','w','Position', [50 100 1800 600]); 
t = tiledlayout(1, 3, 'TileSpacing', 'loose', 'Padding', 'compact');

% Global Colormap for Maps 1 and 2
cmap_maps = [
    190/255 194/255 141/255;   % elimination
    227/255 198/255 143/255;   % equilibrium
    229/255 183/255 166/255;   % escape
];

%% =====================================================================
%% PANEL 1: Lifespan vs k1 (Heatmap)
%% =====================================================================
nexttile;
data1 = load('data_figure5a.mat'); 

[X1, Y1] = meshgrid(data1.lifespan_days, data1.k1_grid);
surf(X1, Y1, zeros(size(data1.Z1)), double(data1.Z1), 'EdgeColor','none');
% contour(X1, Y1, double(data1.Z1), [1.5 2.5], 'k', 'LineWidth', 1.5);
view(2); colormap(gca, cmap_maps); clim([1 3]);
set(gca, 'YScale', 'log', 'XScale', 'linear');
xlim([20, 200]);

xlabel('Memory cell lifespan [days], $1/\epsilon$', 'FontSize', 18);
ylabel('Blast proliferation rate $k_1$ [days$^{-1}]$', 'FontSize', 18);
title('Outcome vs memory cell lifespan', 'FontSize', 20);

hold on;
xline(1/0.014, 'k--', 'LineWidth', 2);
plot([1/0.014, 1/0.014, 1/0.014], [0.06, 0.1, 0.16], 'r*', 'MarkerSize', 10, 'LineWidth', 2);
yline(0.45, 'k-', 'LineWidth', 1.5);
text(25, 0.6, '$k_1^{MRD}$', 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold');

% Regional labels
text(166, 0.008, 'elimination', 'FontSize', 14);
text(166, 0.22, 'elimination', 'FontSize', 14);
text(166, 0.12, 'equilibrium', 'FontSize', 14);
text(176, 0.35, 'escape', 'FontSize', 14);

% --- Add this inside your heatmap plotting section ---
hold on;

% 1. Smoothing & Upsampling (as discussed before)
X_fine = linspace(min(data1.lifespan_days), max(data1.lifespan_days), 400);
Y_fine = logspace(log10(min(data1.k1_grid)), log10(max(data1.k1_grid)), 400);
[XF, YF] = meshgrid(X_fine, Y_fine);

% 2. MASK A: Elimination vs. Everything Else (The Bottom/Outer Boundary)
% Treat Equilibrium (2) and Escape (3) as the same (value = 2)
Z_elim_vs_all = double(data1.Z1);
Z_elim_vs_all(Z_elim_vs_all >= 2) = 2; 
Z_elim_smooth = imgaussfilt(Z_elim_vs_all, 2.0, 'Padding', 'replicate');
ZF_elim = interp2(X1, Y1, Z_elim_smooth, XF, YF, 'cubic');

% Draw the SINGLE line for Elimination boundary
contour(XF, YF, ZF_elim, [1.5 1.5], 'k', 'LineWidth', 2);

% 3. MASK B: Equilibrium vs. Escape (The Inner Boundary)
% We only care about the border between 2 and 3
Z_equil_vs_esc = double(data1.Z1);
Z_equil_smooth = imgaussfilt(Z_equil_vs_esc, 2.0, 'Padding', 'replicate');
ZF_equil = interp2(X1, Y1, Z_equil_smooth, XF, YF, 'cubic');

% Draw the line specifically between 2 and 3
% This will only appear where an Equilibrium zone actually exists
contour(XF, YF, ZF_equil, [2.5 2.5], 'k', 'LineWidth', 2);

hold off;

hold off;

%% =====================================================================
%% PANEL 2: Bhalf vs k1 (Heatmap)
%% =====================================================================
nexttile;
data2 = load('data_figure5b.mat');

[X2, Y2] = meshgrid(data2.Bhalf_grid, data2.k1_grid);
surf(X2, Y2, zeros(size(data2.Z2)), double(data2.Z2), 'EdgeColor','none');
% contour(X2, Y2, double(data2.Z2), [1.5 2.5], 'k', 'LineWidth', 1.5);
view(2); colormap(gca, cmap_maps); clim([1 3]);
set(gca, 'YScale', 'log', 'XScale', 'log');

xlabel('Threshold blast [cell number], $B_{1/2}$', 'FontSize', 18);
ylabel('Blast proliferation rate $k_1$ [days$^{-1}]$', 'FontSize', 18);
title('Outcome vs effector / memory partitionning', 'FontSize', 20);

hold on;
xline(10^9, 'k--', 'LineWidth', 2);
plot(10^10, 0.2, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
yline(0.45, 'k-', 'LineWidth', 1.5);
text(1.3*10^7, 0.6, '$k_1^{MRD}$', 'FontSize', 14, 'Interpreter', 'latex');

% Regional labels 
text(1.8*10^10, 0.24, 'elimination', 'FontSize', 14);
text(1.8*10^10, 0.12, 'equilibrium', 'FontSize', 14);
text(1.2 * 10^7, 0.0025, 'elimination', 'FontSize', 14);
text(3.1 * 10^10, 0.0025, 'escape', 'FontSize', 14);

hold off;

%% =====================================================================
%% PANEL 3: Time Series (PDMP Simulation)
%% =====================================================================
nexttile;
% Simulation Parameters
P_ts = default_params_table1();
P_ts.k1 = 0.18; 
Lambda  = 1e2;
MRD = 1e6; ymin = 1; ymax = 1e12;

% Run the exact simulation code provided
sim_ts = simulate_beam_exact_pdmp(P_ts, 1000, 0.1, Lambda);

% Call the plot helper
plot_timeseries_panel(gca, sim_ts, MRD, ymin, ymax, 1);
title('Elimination after relapse', 'FontSize', 20);

%% =====================================================================
%% GLOBAL STYLING & EXPORT
%% =====================================================================
% Apply LaTeX and Font settings to all axes in the figure
all_axes = findall(gcf, 'type', 'axes');
for ax = all_axes'
    set(ax, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 14, 'Box', 'on', 'Layer', 'top', 'GridAlpha', 0.3);
    set(ax.Title, 'Interpreter', 'latex');
    set(ax.XLabel, 'Interpreter', 'latex');
    set(ax.YLabel, 'Interpreter', 'latex');
end

% Set default interpreters for any text called later
set(gcf, 'DefaultTextInterpreter', 'latex');

% 1. Set units to physical measurements
set(gcf, 'Units', 'inches');

% 2. Define [left, bottom, width, height]
% This makes the figure 15 inches wide and 5 inches tall
set(gcf, 'Position', [1, 1, 20, 5]); 

% 3. Crucial: Ensure the PaperPosition matches for older export methods
% (exportgraphics usually handles this, but it's good practice)
set(gcf, 'PaperPositionMode', 'auto');

% % Add a shared colorbar for the first two maps
% cb = colorbar(all_axes(2)); % Attach to the second plot
% cb.Layout.Tile = 'east'; 
% cb.Ticks = [1 2 3];
% cb.TickLabels = {'elimination','equilibrium','escape'};
% cb.TickLabelInterpreter = 'latex';

% Exporting
exportgraphics(gcf, 'Combined_Figure.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%% =====================================================================
%% REQUIRED FUNCTIONS (Place these at the end of your .m file)
%% =====================================================================

function P = default_params_table1()
    P.k1 = 0.2; P.k2 = 2.9e-10; P.k3 = 1e-9; P.k4 = 0.1;
    P.gamma = 0.3; P.delta = 0.2; P.eps = 0.014; P.K = 1e12;
    P.N = 6; P.Bhalf = 1e10; P.CAR0 = 4.1e8; P.M0 = .468 * P.CAR0;
    P.A0 = 2000; P.E0each = (.532 * P.CAR0)/P.N; P.B0 = 2e11;
end

% ... Include simulate_beam_exact_pdmp, hybrid_step_exact_pdmp, 
% ... pdmp_ode_exact, pdmp_event, get_discrete_propensities, 
% ... beam_ode_full, and plot_timeseries_panel from your original code here.


% =====================================================================
% Hybrid simulation main loop
% =====================================================================
function out = simulate_beam_exact_pdmp(P, Tend, dt, Lambda)

N  = P.N;
tgrid = (0:dt:Tend).';
ng = numel(tgrid);

Blog = zeros(ng,1); Elog = zeros(ng,1); Mlog = zeros(ng,1); Tlog = zeros(ng,1);

B = P.B0; E = ones(N,1)*P.E0each; A = P.A0; M = P.M0;

k = 1;
Blog(k) = B; Elog(k) = sum(E); Mlog(k) = M + A; Tlog(k) = sum(E) + M + A;

t = 0;
for k = 2:ng
    tNext = tgrid(k);
    h = tNext - t; 

    if B >= Lambda && A >= Lambda && M >= Lambda && all(E >= Lambda)
        % Deterministic step
        y0 = [B; E; A; M];
        opts = odeset('RelTol',1e-7,'AbsTol',1e-10);
        [~,yy] = ode15s(@(tt,yy) beam_ode_full(tt,yy,P), [t tNext], y0, opts);
        yend = yy(end,:).';
        B = max(yend(1),0); E = max(yend(2:1+N),0); A = max(yend(2+N),0); M = max(yend(3+N),0);
    else
        % Exact PDMP Step
        [B,E,A,M] = hybrid_step_exact_pdmp(P, B, E, A, M, h, Lambda);
    end

    t = tNext;
    % disp(t)
    Blog(k) = B; Elog(k) = sum(E); Mlog(k) = M + A; Tlog(k) = sum(E) + M + A;

    if B <= 0, B = 0; end
end

out.t = tgrid; out.B = Blog; out.E = Elog; out.M = Mlog; out.T = Tlog;      
out.eliminated = any(Blog <= 0.5); 
end

% =====================================================================
% EXACT PDMP Hybrid Step (Integral Method)
% =====================================================================
function [B,E,A,M] = hybrid_step_exact_pdmp(P, B, E, A, M, h, Lambda)

N = P.N;
tloc = 0;

% Determine which states are discrete for this window and round them
is_B_disc = (B < Lambda); if is_B_disc, B = round(B); end
is_A_disc = (A < Lambda); if is_A_disc, A = round(A); end
is_M_disc = (M < Lambda); if is_M_disc, M = round(M); end
is_E_disc = (E < Lambda);  
for i=1:N
    if is_E_disc(i), E(i) = round(E(i)); end
end

% Reaction Flags (1 = Continuous, 0 = Discrete)
f_flags = zeros(6+N, 1);
f_flags(1) = ~is_B_disc;                        
f_flags(2) = ~is_B_disc;                        
f_flags(3) = ~(is_M_disc || is_A_disc);         
f_flags(4) = ~(is_A_disc || is_M_disc);         
f_flags(5) = ~(is_A_disc || is_E_disc(1));      
f_flags(6) = ~is_M_disc;                        
for i=1:N-1
    f_flags(6+i) = ~(is_E_disc(i) | is_E_disc(i+1)); 
end
f_flags(6+N) = ~is_E_disc(N);                   

while tloc < h
    % 1. Draw the exact integral threshold R = -log(u)
    R = -log(rand);
    
    % 2. Setup ODE with event locator
    % State vector: y = [B; E(1..N); A; M; I] where I is the accumulated propensity integral
    y0 = [B; E; A; M; 0];
    opts = odeset('RelTol',1e-7,'AbsTol',1e-10, 'Events', @(t,y) pdmp_event(t,y,R,P,f_flags));
    
    % 3. Integrate until event fires (I == R) or we reach end of step (h - tloc)
    [tt, yy, te, ye, ie] = ode15s(@(t,y) pdmp_ode_exact(t,y,P,f_flags), [0, h - tloc], y0, opts);
    
    % Update populations to the end of the integration period
    yend = yy(end,:).';
    B = max(yend(1),0); E = max(yend(2:1+N),0); A = max(yend(2+N),0); M = max(yend(3+N),0);
    
    if ~isempty(te)
        % 4. A discrete event fired! Update local time.
        tloc = tloc + te(1);
        
        % 5. Compute instantaneous conditional probabilities P[j = i | T = t+tau]
        props = get_discrete_propensities([B; E; A; M], P, f_flags);
        a0 = sum(props);
        
        if a0 > 0
            u = rand * a0;
            cumsum_props = cumsum(props);
            idx = find(cumsum_props > u, 1);
            
            % 6. Apply exact integer stoichiometry
            switch idx
                case 1, B = B + 1;
                case 2, B = max(B - 1, 0);
                case 3, M = max(M - 1, 0); A = A + 1;
                case 4, A = max(A - 1, 0); M = M + 2;
                case 5, A = max(A - 1, 0); E(1) = E(1) + 2;
                case 6, M = max(M - 1, 0);
                case 6 + N, E(N) = max(E(N) - 1, 0);
                otherwise 
                    i = idx - 6; 
                    E(i) = max(E(i) - 1, 0); E(i+1) = E(i+1) + 2;
            end
        end
    else
        % Reached the end of the time step without any discrete events firing
        tloc = h;
    end
end
end

% =====================================================================
% ODE: Exact Subsystem with Propensity Integral
% =====================================================================
function dydt = pdmp_ode_exact(~, y, P, f_flags)
N = P.N;
B = y(1); E = y(2:1+N); A = y(2+N); M = y(3+N);

Etot = sum(E); frac = B/(P.Bhalf + B);

% Continuous population changes (flags = 1)
dB = f_flags(1) * P.k1*B*(1 - B/P.K) - f_flags(2) * P.k2*B*Etot;

dE = zeros(N,1);
dE(1) = f_flags(5) * 2*P.k4*A*frac - f_flags(7) * P.gamma*E(1);
for i = 2:N-1
    dE(i) = f_flags(6+i-1) * 2*P.gamma*E(i-1) - f_flags(6+i) * P.gamma*E(i);
end
dE(N) = f_flags(6+N-1) * 2*P.gamma*E(N-1) - f_flags(6+N) * P.delta*E(N);

dA = f_flags(3) * P.k3*M*B - f_flags(4) * P.k4*A*(1-frac) - f_flags(5) * P.k4*A*frac;
dM = -f_flags(3) * P.k3*M*B + f_flags(4) * 2*P.k4*A*(1-frac) - f_flags(6) * P.eps*M;

% The Integral of discrete propensities (flags = 0)
props = get_discrete_propensities([B; E; A; M], P, f_flags);
dI = sum(props);

dydt = [dB; dE; dA; dM; dI];
end

% =====================================================================
% Event function to halt ODE exactly when Integral == R
% =====================================================================
function [value, isterminal, direction] = pdmp_event(~, y, R, ~, ~)
I = y(end);
value = I - R;      % Triggers when I crosses R
isterminal = 1;     % Halt integration
direction = 1;      % Only trigger when crossing from below
end

% =====================================================================
% Helper: Calculate Propensities for Discrete Reactions
% =====================================================================
function props = get_discrete_propensities(y, P, f_flags)
N = P.N;
B = y(1); E = y(2:1+N); A = y(2+N); M = y(3+N);

Etot = sum(E); frac = B/(P.Bhalf + B);
props = zeros(6+N, 1);

if ~f_flags(1), props(1) = max(P.k1 * B * (1 - B/P.K), 0); end
if ~f_flags(2), props(2) = max(P.k2 * B * Etot, 0); end
if ~f_flags(3), props(3) = max(P.k3 * M * B, 0); end
if ~f_flags(4), props(4) = max(P.k4 * A * (1 - frac), 0); end
if ~f_flags(5), props(5) = max(P.k4 * A * frac, 0); end
if ~f_flags(6), props(6) = max(P.eps * M, 0); end

for i = 1:N-1
    if ~f_flags(6+i)
        props(6+i) = max(P.gamma * E(i), 0);
    end
end
if ~f_flags(6+N), props(6+N) = max(P.delta * E(N), 0); end

end

% =====================================================================
% ODE: Full System (Deterministic fallback)
% =====================================================================
function dydt = beam_ode_full(~, y, P)
N = P.N;
B = y(1); E = y(2:1+N); A = y(2+N); M = y(3+N);
Etot = sum(E); frac = B/(P.Bhalf + B);

dB = P.k1*B*(1 - B/P.K) - P.k2*B*Etot;
dE = zeros(N,1);
dE(1) = 2*P.k4*A*frac - P.gamma*E(1);
for i = 2:N-1
    dE(i) = P.gamma*(2*E(i-1) - E(i));
end
dE(N) = 2*P.gamma*E(N-1) - P.delta*E(N);
dA = P.k3*M*B - P.k4*A*(1-frac) - P.k4*A*frac;
dM = -P.k3*M*B + 2*P.k4*A*(1-frac) - P.eps*M;

dydt = [dB; dE; dA; dM];
end

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

set(ax, 'YScale','log'); xlim(ax, [0 600]); ylim(ax, [ymin ymax]);
xlabel(ax, 'Time [days]'); ylabel(ax, 'Cell Number');
grid(ax,'on'); set(ax,'Box','on');
if position==1
    text(500, 2.5e6, 'MRD positive', 'FontSize', 9, 'Color',[0.3 0.3 0.3]);
    text(500, 6e5,  'MRD negative', 'FontSize', 9, 'Color',[0.3 0.3 0.3]);
else
    text(500, 2.5e6, 'MRD positive', 'FontSize', 9, 'Color',[0.3 0.3 0.3]);
    text(500, 6e9,  'MRD negative', 'FontSize', 9, 'Color',[0.3 0.3 0.3]);
end
legend({'Blast','Effector','Memory','Total'}, 'Location','northeast', 'Box','on', 'NumColumns', 2);
hold off;
end