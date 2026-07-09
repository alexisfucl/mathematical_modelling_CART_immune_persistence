
plot_figure7_from_mat('figure7_oldTN_tauLeapFlux_ACD_20260515_101629.mat', true)
% fig = plot_figure7_from_mat('yourfile.mat', false);

function fig = plot_figure7_from_mat(matFile, savePng)
%PLOT_FIGURE7_FROM_MAT  Recreate Figure 7 (a,c,d) plots from a saved .mat file.
% Uses a 2-color gradient:
%   p=0 (escape)      -> Fig5 escape red/pink
%   p=1 (elimination) -> Fig5 elimination green
% And applies LaTeX + larger styling everywhere.

if nargin < 2, savePng = false; end
if ~isfile(matFile)
    error('MAT file not found: %s', matFile);
end

S = load(matFile);

req = {'P_elim_a','alpha_grid','tau_a', ...
       'P_elim_c','gamma_grid','tau_c', ...
       'P_elim_d','delta_grid','tau_d'};
for k = 1:numel(req)
    if ~isfield(S, req{k})
        error('Missing field "%s" in MAT file: %s', req{k}, matFile);
    end
end

smoothSigma = 1.2;
if isfield(S,'smoothSigma'), smoothSigma = S.smoothSigma; end

% --- 2-color gradient using Fig 5 endpoints ---
cmap_prob = cmap_fig5_escape_to_elim(256);

% ---- Global styling (LaTeX + font sizes) ----
FS_tick  = 16;
FS_label = 20;
FS_title = 22;
FS_cb    = 18;

fig = figure('Color','w','Position',[60 120 1700 560]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% -------- (a)
ax1 = nexttile; hold(ax1,'on');
cb1 = plot_probmap_with_half_contour(ax1, S.tau_a, S.alpha_grid, S.P_elim_a, smoothSigma, cmap_prob);
set(ax1,'YScale','log');
xlabel(ax1,'Mean time in the niche $1/m$ [days]','FontSize',FS_label,'Interpreter','latex');
ylabel(ax1,'Entry rate $\alpha$ [days$^{-1}$]','FontSize',FS_label,'Interpreter','latex');
title(ax1,'(a)','FontSize',FS_title,'Interpreter','latex');
style_axes(ax1, FS_tick);

cb1.Label.String = 'Elimination probability';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize = FS_cb;
cb1.TickLabelInterpreter = 'latex';
cb1.FontSize = FS_tick;
cb1.Ticks = [0 0.5 1];
hold(ax1,'off');

% -------- (c)
ax2 = nexttile; hold(ax2,'on');
cb2 = plot_probmap_with_half_contour(ax2, S.tau_c, S.gamma_grid, S.P_elim_c, smoothSigma, cmap_prob);
set(ax2,'YScale','log');
xlabel(ax2,'Mean time in the niche $1/m$ [days]','FontSize',FS_label,'Interpreter','latex');
ylabel(ax2,'Effector division rate $\gamma$ [days$^{-1}$]','FontSize',FS_label,'Interpreter','latex');
title(ax2,'(b)','FontSize',FS_title,'Interpreter','latex');
style_axes(ax2, FS_tick);

cb2.Label.String = 'Elimination probability';
cb2.Label.Interpreter = 'latex';
cb2.Label.FontSize = FS_cb;
cb2.TickLabelInterpreter = 'latex';
cb2.FontSize = FS_tick;
cb2.Ticks = [0 0.5 1];
hold(ax2,'off');

% -------- (d)
ax3 = nexttile; hold(ax3,'on');
cb3 = plot_probmap_with_half_contour(ax3, S.tau_d, S.delta_grid, S.P_elim_d, smoothSigma, cmap_prob);
set(ax3,'YScale','log');
xlabel(ax3,'Mean time in the niche $1/m$ [days]','FontSize',FS_label,'Interpreter','latex');
ylabel(ax3,'Effector death rate $\delta$ [days$^{-1}$]','FontSize',FS_label,'Interpreter','latex');
title(ax3,'(c)','FontSize',FS_title,'Interpreter','latex');
style_axes(ax3, FS_tick);

cb3.Label.String = 'Elimination probability';
cb3.Label.Interpreter = 'latex';
cb3.Label.FontSize = FS_cb;
cb3.TickLabelInterpreter = 'latex';
cb3.FontSize = FS_tick;
cb3.Ticks = [0 0.5 1];
hold(ax3,'off');

% Apply default LaTeX text interpreter for any later text calls
set(fig,'DefaultTextInterpreter','latex');

% Save PNG if requested
if savePng
    [p,n,~] = fileparts(matFile);
    outPng = fullfile(p, [n '_replot.png']);
    print(fig, outPng, '-dpng', '-r300');
    fprintf('Saved PNG: %s\n', outPng);
end

end

%% ========================= helpers =========================
function cb = plot_probmap_with_half_contour(ax, xVals, yVals, P, sigmaPix, cmap)

[X,Y] = meshgrid(xVals, yVals);
surf(ax, X, Y, zeros(size(P)), P, 'EdgeColor','none');
view(ax, 2);

colormap(ax, cmap);
caxis(ax, [0 1]);
grid(ax,'on'); box(ax,'on');
set(ax,'Layer','top');

Ps = gaussian_smooth_2d(P, sigmaPix);
contour(ax, xVals, yVals, Ps, [0.5 0.5], 'k-', 'LineWidth', 2.5);

cb = colorbar(ax);
end

function Ps = gaussian_smooth_2d(P, sigma)
if sigma <= 0
    Ps = P; return;
end
rad = max(1, ceil(3*sigma));
x = (-rad:rad);
g = exp(-(x.^2)/(2*sigma^2));
g = g / sum(g);
Ps = conv2(conv2(P, g, 'same'), g', 'same');
end

function cmap = cmap_fig5_escape_to_elim(n)
col_elim = [0.70 0.78 0.55];  % elimination green
col_esc  = [0.88 0.74 0.70];  % escape red/pink
t = linspace(0,1,n).';
cmap = (1-t).*col_esc + t.*col_elim;
cmap = max(min(cmap,1),0);
end

function style_axes(ax, FS_tick)
set(ax, ...
    'TickLabelInterpreter','latex', ...
    'FontName','Times', ...
    'FontSize',FS_tick, ...
    'Box','on', ...
    'Layer','top', ...
    'GridAlpha',0.25);
end