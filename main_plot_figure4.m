function fig3_beam_reproduce()
%FIG3_BEAM_REPRODUCE Reproduce Figure 3 of the BEAM paper
% Implements the mathematically exact Integral Method for the PDMP hybrid
% model. Continuous populations dictate time-varying discrete propensities.

rng(10,'twister');

% =========================
% Parameters (Table 1 baseline)
% =========================
P = default_params_table1();

Tend    = 1000;     
dt      = 0.1;      
Lambda  = 1e2;     
MRD     = 1e6;     
ymin    = 1e0;
ymax    = 1e12;

% k1_vals   = [0.06, 0.10, 0.16]; % values for k2 = 2.9
k1_vals   = [0.04, 0.10, 0.157];

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

% This forces the tiles to fill the entire rectangular area 
% defined by the figure size.
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
% 1. Set the figure to be wide (e.g., 10 inches wide, 4 inches tall)
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1, 1, 20, 4]); % [left, bottom, width, height]

% 2. Ensure the axes expand to fill that new rectangular shape
% If using subplot:
set(findall(gcf,'type','axes'), 'PlotBoxAspectRatioMode', 'auto');

% 3. The magic command for LaTeX users
% This crops the PDF to the actual bounding box of the content
exportgraphics(gcf, 'Figure4.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
end

% =====================================================================
% Defaults
% =====================================================================
function P = default_params_table1()
P.k1     = 0.2;
P.k2     = 2.6e-10;     
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
P.A0     = 2000; % Start at threshold to bypass massive initial flux                    
P.E0each = (.532 * P.CAR0)/P.N; 
P.B0     = 2e11;
end

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