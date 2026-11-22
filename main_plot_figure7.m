plot_figure7_from_mat('figure7_data.mat', true);

function fig = plot_figure7_from_mat(matFile, savePng)
if nargin < 2, savePng = false; end
if ~isfile(matFile)
    error('MAT file not found: %s', matFile);
end

S = load(matFile);

req = {'tau','alpha_grid','gamma_grid','delta_grid','P_elim_a','P_elim_b','P_elim_c'};
for k = 1:numel(req)
    if ~isfield(S, req{k})
        error('Missing field "%s" in MAT file: %s', req{k}, matFile);
    end
end

smoothSigma = 0;   % the 0.5 contour is drawn on the raw probabilities

cmap_prob = cmap_fig5_escape_to_elim(256);


FS_tick  = 16;
FS_label = 20;
FS_title = 22;
FS_cb    = 18;

tauMax = 50; gammaMin = 0.006; deltaMax = 0.2;

fig = figure('Color','w','Position',[60 120 1700 460]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% -------- (a)
ax1 = nexttile; hold(ax1,'on');
cb1 = plot_probmap_with_half_contour(ax1, S.tau, S.alpha_grid, S.P_elim_a, smoothSigma, cmap_prob);
set(ax1,'YScale','log'); xlim(ax1, [0 tauMax]); ylim(ax1, [1e-13 1e-9]);
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

% -------- (b)
ax2 = nexttile; hold(ax2,'on');
cb2 = plot_probmap_with_half_contour(ax2, S.tau, S.gamma_grid, S.P_elim_b, smoothSigma, cmap_prob);
set(ax2,'YScale','log'); xlim(ax2, [0 tauMax]); ylim(ax2, [gammaMin S.gamma_grid(end)]);
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

% -------- (c)
ax3 = nexttile; hold(ax3,'on');
cb3 = plot_probmap_with_half_contour(ax3, S.tau, S.delta_grid, S.P_elim_c, smoothSigma, cmap_prob);
set(ax3,'YScale','log'); xlim(ax3, [0 tauMax]); ylim(ax3, [S.delta_grid(1) deltaMax]);
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

set(fig,'DefaultTextInterpreter','latex');
if savePng
    for a = [ax1 ax2 ax3], a.Toolbar.Visible = 'off'; end   
    drawnow;
    exportgraphics(fig, 'Figure_7.pdf', 'ContentType', 'vector');
    exportgraphics(fig, 'Figure_7.png', 'Resolution', 300);
end
end


function cb = plot_probmap_with_half_contour(ax, xVals, yVals, P, sigmaPix, cmap)
[X,Y] = meshgrid(cell_edges(xVals, false), cell_edges(yVals, true));
surf(ax, X, Y, zeros(size(X)), P([1:end end], [1:end end]), 'EdgeColor','flat', 'LineWidth',2);
view(ax, 2);
xlim(ax, [0 xVals(end)]); ylim(ax, [yVals(1) yVals(end)]);

colormap(ax, cmap);
caxis(ax, [0 1]);
grid(ax,'on'); box(ax,'on');
set(ax,'Layer','top');

Ps = gaussian_smooth_2d(P, sigmaPix);
contour(ax, xVals, yVals, Ps, [0.5 0.5], 'k-', 'LineWidth', 2.5);

cb = colorbar(ax);
end

function e = cell_edges(v, logScale)
v = v(:)';
if logScale, v = log10(v); end
h = diff(v);
e = [v(1) - h(1)/2, v(1:end-1) + h/2, v(end) + h(end)/2];
if logScale, e = 10.^e; end
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
col_elim = [0.70 0.78 0.55];  
col_esc  = [0.88 0.74 0.70];  
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
