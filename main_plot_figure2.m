% plot_figure2_bifurcation.m
% Reproduce Figure 2 bifurcation diagram (steady states b vs epsilon)
% Direct plotting of equilibrium branches (no ODE simulation).

clear; close all; clc;

%% Global LaTeX interpreters (axes/legend/text)
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Parameters (choose to match the visual in the screenshot)
b_half = 0.7;   % b_{1/2} (matches the red-branch intercept near 0.5)
k3     = 5.0;   % scale so epsilon spans [0,1] in the figure

eps_max  = 1.0;

%% Fold location: epsilon* = (b_{1/2} k3) / (3 + 2 sqrt(2))
eps_fold = (b_half*k3) / (3 + 2*sqrt(2));

%% Build epsilon grids
% - For constant branches (b=0, b=1): full [0,1]
eps_all = linspace(0, eps_max, 1500);

% - For bifurcating branches: only defined for eps in [0, eps_fold]
%   Add strong refinement near eps_fold to avoid "blank space".
Ncoarse = 1200;
Nref    = 2500;

eps_branch_coarse = linspace(0, eps_fold, Ncoarse);

% cluster points extremely close to the fold from the left:
% eps = eps_fold - delta, with delta spanning many decades
delta = eps_fold * logspace(-10, -2, Nref);
eps_branch_refined = eps_fold - delta;
eps_branch_refined(eps_branch_refined < 0) = [];

% include the fold point exactly
eps_branch = unique([eps_branch_coarse, eps_branch_refined, eps_fold]);
eps_branch = sort(eps_branch);

%% Compute b*_± from Eq. (2.8a)
D = (b_half*k3 - eps_branch).^2 - 4*b_half*k3.*eps_branch;  % discriminant
% numerical guard (tiny negatives from roundoff)
D(D < 0 & D > -1e-14) = 0;

b_plus  = ((b_half*k3 - eps_branch) + sqrt(D)) ./ (2*k3);  % unstable (red)
b_minus = ((b_half*k3 - eps_branch) - sqrt(D)) ./ (2*k3);  % stable   (blue)

%% Constant steady states
b0 = zeros(size(eps_all));  % b0* = 0 (unstable)
b1 = ones(size(eps_all));   % b1* = 1 (stable)

%% Plot
figure('Color','w'); hold on;

% Stable (dark blue)
darkBlue = [0 0.2 0.6];
hStable1 = plot(eps_all, b1, '-', 'Color', darkBlue, 'LineWidth', 2);
hStable2 = plot(eps_branch, b_minus, '-', 'Color', darkBlue, 'LineWidth', 2);

% Unstable (dark red)
darkRed = [.75 0 0];
hUnst1 = plot(eps_all, b0, '-', 'Color', darkRed, 'LineWidth', 2);
hUnst2 = plot(eps_branch, b_plus, '-', 'Color', darkRed, 'LineWidth', 2);
% Fold marker (vertical dashed line)
plot([eps_fold eps_fold], [0, max(b_plus)*1.05], 'k--', 'LineWidth', 1.2);

%% Axes formatting to match screenshot
xlim([0 1]);
ylim([-0.05 1.15]);

% Only show 0 and 1 ticks like the screenshot
xticks([0 1]);
yticks([0 1]);
yticklabels({'$0$' ,'$K$'});

xlabel('$\epsilon$', 'FontSize', 18);
ylabel('$B$', 'FontSize', 22, 'Rotation', 0);  % horizontal "facing down" look
ax = gca;
ax.YLabel.HorizontalAlignment = 'right';
ax.YLabel.VerticalAlignment   = 'middle';

ax.FontSize  = 14;
ax.LineWidth = 1;
ax.Box       = 'on';
ax.TickDir   = 'out';   % clean ticks like typical publication plots
ax.TickLength = [0.02 0.02];

%% Legend 
% single-line legend: place entries side-by-side
lg = legend([hStable1, hUnst1], {'Stable','Unstable'}, 'Location','north', ...
    'Orientation','horizontal', 'Box','on');
lg.Box = 'on';

%% Annotations (LaTeX)
text(0.92, 1.05, '$B_1^\ast$', 'FontSize', 14);
text(0.92, 0.05, '$B_0^\ast$', 'FontSize', 14);

text(-0.08, b_half, '$B_{1/2}$', 'FontSize', 14);

% Put b*_{+} and b*_{-} labels near the branches
eps_lab = 0.55*eps_fold;
[~, idx] = min(abs(eps_branch - eps_lab));
text(eps_branch(idx)+0.02, b_plus(idx)+0.03,  '$B_{+}^\ast$', 'FontSize', 14);
text(eps_branch(idx)+0.02, b_minus(idx)+0.07, '$B_{-}^\ast$', 'FontSize', 14);

% Fold formula under dashed line (as in screenshot)
txtFold = '$\epsilon=\frac{B_{1/2}k_3}{3+2\sqrt{2}}$';
text(eps_fold-0.1, -0.12, txtFold, 'FontSize', 10);

hold off;

