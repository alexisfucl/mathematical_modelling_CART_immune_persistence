% plot_figure2_bifurcation.m
% Bifurcation diagram (steady states B vs epsilon)
% Direct plotting of equilibrium branches (no ODE simulation).

clear; close all; clc;

%% Parameters (dimensional)
p.K     = 1e12;                 % cells
p.H     = 0.7*p.K;              % B_{1/2}, cells: enlarged for visibility
p.k1    = 0.05;                 % day^-1: illustrative stable regime
p.k2    = 4e-10;                % cell^-1 day^-1
p.k3    = 5/p.K;                % cell^-1 day^-1: k3*K = 5 day^-1
p.k4    = 0.1;                  % day^-1
p.gamma = 0.3;                  % day^-1
p.delta = 0.2;                  % day^-1
p.N     = 6;
eps_max = 1.0;                  % day^-1

eps_fold = p.k3*p.H/(3 + 2*sqrt(2));
B_fold   = p.H*(sqrt(2) - 1);
assert(p.H < p.K && eps_max > eps_fold);

%% Equilibrium branches: refine near both ends, excluding the fold itself
fraction = unique([logspace(-6,-2,251), linspace(0.01,0.99,5001), ...
                   1-logspace(-10,-2,251)]);
eps_branch = eps_fold*fraction;
n = numel(eps_branch);
B_minus = zeros(1,n); B_plus = zeros(1,n);
lambda_minus = zeros(1,n); lambda_plus = zeros(1,n);

for j = 1:n
    [xm,Jm] = coexistence(eps_branch(j),-1,p);
    [xp,Jp] = coexistence(eps_branch(j),+1,p);
    B_minus(j) = xm(1);
    B_plus(j)  = xp(1);
    lambda_minus(j) = max(real(eig(Jm)));
    lambda_plus(j)  = max(real(eig(Jp)));
end

% Stop rather than mislabel the branches if the parameters are changed.
assert(all(lambda_minus < 0), ...
    'The lower branch is not stable over the whole sampled interval.');
assert(all(lambda_plus > 0), ...
    'The upper branch is not unstable over the whole sampled interval.');

% The fold is nonhyperbolic: one zero eigenvalue, all others negative.
[~,Jfold] = coexistence(eps_fold,-1,p);
fold_eigenvalues = eig(Jfold);
[~,izero] = min(abs(fold_eigenvalues));
other_eigenvalues = fold_eigenvalues;
other_eigenvalues(izero) = [];
assert(abs(fold_eigenvalues(izero)) < 1e-9);
assert(all(real(other_eigenvalues) < 0));

%% Plot: raw B values, so no division by B_{1/2} and no broken axis
blue = [0 0.2 0.6];
red  = [0.75 0 0];
fig = figure('Color','w','Units','inches','Position',[1 1 5.3 4.9]);
ax = axes(fig,'Position',[0.20 0.17 0.75 0.72]);
ax.Toolbar.Visible = 'off';
disableDefaultInteractivity(ax);
hold(ax,'on');

hStable = plot(ax,[0 eps_max],[p.K p.K],'-', ...
    'Color',blue,'LineWidth',2);
plot(ax,[eps_branch eps_fold],[B_minus B_fold],'-', ...
    'Color',blue,'LineWidth',2);
hUnstable = plot(ax,[0 eps_max],[0 0],'-', ...
    'Color',red,'LineWidth',2);
plot(ax,[eps_branch eps_fold],[B_plus B_fold],'-', ...
    'Color',red,'LineWidth',2);

plot(ax,[eps_fold eps_fold],[0 B_fold],'k:','LineWidth',2);
plot(ax,eps_fold,B_fold,'o','MarkerSize',6, ...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.1);

xlim(ax,[0 eps_max]);
ylim(ax,[-0.05 1.25]*p.K);
xticks(ax,[0 eps_fold eps_max]);
xticklabels(ax,{'$0$','','$1$'});
% Separate text allows a larger fold label without enlarging every tick.
text(ax,eps_fold/eps_max,-0.037,'$\epsilon_{\mathrm{f}}$', ...
    'Units','normalized','Interpreter','latex','FontSize',17, ...
    'HorizontalAlignment','center','VerticalAlignment','top', ...
    'Clipping','off');
yticks(ax,[0 p.H p.K]);
yticklabels(ax,{'$0$','$B_{1/2}$','$K$'});
ax.YAxis.Exponent = 0;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 13;
ax.LineWidth = 1;
ax.Box = 'on';
ax.TickDir = 'out';
ax.TickLength = [0.018 0.018];
xlabel(ax,'Memory CAR T-cell death rate, $\epsilon$ (day$^{-1}$)', ...
    'Interpreter','latex','FontSize',13);
ylabel(ax,'Equilibrium blast burden, $B^*$', ...
    'Interpreter','latex','FontSize',14,'Rotation',90);

legend(ax,[hStable hUnstable],{'Stable','Unstable'}, ...
    'Interpreter','latex','Location','northoutside', ...
    'Orientation','horizontal','Box','off','FontSize',13);
text(ax,0.93,1.025*p.K,'$B_1^*$','Interpreter','latex', ...
    'FontSize',14,'VerticalAlignment','bottom');
text(ax,0.93,0.025*p.K,'$B_0^*$','Interpreter','latex', ...
    'FontSize',14,'VerticalAlignment','bottom');
[~,jlabel] = min(abs(eps_branch-0.48*eps_fold));
text(ax,eps_branch(jlabel),B_plus(jlabel)+0.025*p.K, ...
    '$B_+^*$','Interpreter','latex','FontSize',14,'VerticalAlignment','bottom');
text(ax,eps_branch(jlabel),B_minus(jlabel)+0.045*p.K, ...
    '$B_-^*$','Interpreter','latex','FontSize',14,'VerticalAlignment','bottom');
hold(ax,'off');

%% Save the actual MATLAB output and the numerical stability check
outdir = fileparts(mfilename('fullpath'));
if isempty(outdir), outdir = pwd; end
drawnow;
exportgraphics(fig,fullfile(outdir,'Figure2_large_Bhalf.pdf'), ...
    'ContentType','vector','BackgroundColor','white','Padding',12);
exportgraphics(fig,fullfile(outdir,'Figure2_large_Bhalf.png'), ...
    'Resolution',300,'BackgroundColor','white','Padding',12);
savefig(fig,fullfile(outdir,'Figure2_large_Bhalf.fig'));

T = table(eps_branch(:),B_minus(:),B_plus(:), ...
          lambda_minus(:),lambda_plus(:),'VariableNames', ...
          {'epsilon_per_day','B_minus_cells','B_plus_cells', ...
           'max_real_lambda_minus','max_real_lambda_plus'});
writetable(T,fullfile(outdir,'figure2_stability_data.csv'));
save(fullfile(outdir,'figure2_parameters.mat'), ...
    'p','eps_fold','B_fold','fold_eigenvalues');

report = sprintf(['MATLAB %s\n' ...
    'k1 = %.3g day^-1\nB_half = %.3g cells; K = %.3g cells\n' ...
    'k3 = %.6g cell^-1 day^-1\n' ...
    'Fold: epsilon = %.12g day^-1; B = %.12g cells\n' ...
    'Checked %d epsilon values strictly between zero and the fold.\n' ...
    'Lower branch: stable at %d/%d sampled values.\n' ...
    'Upper branch: unstable at %d/%d sampled values.\n' ...
    'Fold: one zero eigenvalue; all other eigenvalues have negative real parts.\n'], ...
    version,p.k1,p.H,p.K,p.k3,eps_fold,B_fold,n, ...
    sum(lambda_minus<0),n,sum(lambda_plus>0),n);
fprintf('%s',report);
fid = fopen(fullfile(outdir,'stability_check.txt'),'w');
assert(fid>=0,'Could not write the stability report.');
fprintf(fid,'%s',report);
fclose(fid);

%% Full coexistence state and Jacobian
function [x,Jscaled] = coexistence(epsilon,branch,p)
    z = p.k3*p.H;
    eps_fold = z/(3+2*sqrt(2));
    D = max(0,(z-epsilon)^2-4*z*epsilon);
    if epsilon == eps_fold
        B = p.H*(sqrt(2)-1);
    elseif branch < 0
        % Equivalent quadratic root, avoiding cancellation near epsilon=0.
        B = 2*p.H*epsilon/(z-epsilon+sqrt(D));
    else
        B = (z-epsilon+sqrt(D))/(2*p.k3);
    end
    N = p.N;
    Q = 2^(N-1)-1+(p.gamma/p.delta)*2^(N-1);
    E1 = (p.k1/p.k2)*(1-B/p.K)/Q;
    E = [E1*2.^(0:N-2), (p.gamma/p.delta)*2^(N-1)*E1];
    A = p.gamma*E1*(p.H+B)/(2*p.k4*B);
    M = p.k4*A/(p.k3*B);
    x = [B;E(:);A;M];

    % State ordering: B, E1,...,EN, A, M.
    ia = N+2; im = N+3;
    J = zeros(N+3);
    J(1,1) = p.k1*(1-2*B/p.K)-p.k2*sum(E);
    J(1,2:N+1) = -p.k2*B;
    J(2,1) = 2*p.k4*A*p.H/(p.H+B)^2;
    J(2,2) = -p.gamma;
    J(2,ia) = 2*p.k4*B/(p.H+B);
    for i = 3:N
        J(i,i-1) = 2*p.gamma;
        J(i,i) = -p.gamma;
    end
    J(N+1,N) = 2*p.gamma;
    J(N+1,N+1) = -p.delta;
    J(ia,1) = p.k3*M;
    J(ia,ia) = -p.k4;
    J(ia,im) = p.k3*B;
    J(im,1) = -p.k3*M-2*p.k4*A*p.H/(p.H+B)^2;
    J(im,ia) = 2*p.k4*p.H/(p.H+B);
    J(im,im) = -(p.k3*B+epsilon);
    % Diagonal similarity transformation preserves all eigenvalues.
    Jscaled = (J.*x.')./x;
end
