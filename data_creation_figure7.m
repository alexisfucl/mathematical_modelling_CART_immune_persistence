function compute_figure7_data_hybrid(varargin)

opt = struct('nReps', 500, 'Tend', 180, 'every', 1, 'chunk', [1 1], ...
             'runDir', 'figure7_hybrid_runs', 'outFile', 'figure7_data_hybrid.mat');
for q = 1:2:numel(varargin), opt.(varargin{q}) = varargin{q+1}; end

P = default_params_table1();
dt = 0.1;  Lambda = 1e3;                          

tau        = linspace(0.01, 50, 20);              
alpha_grid = 10.^(-9 - (26:-1:0)*3/19);          
gamma_grid = logspace(-3, -0.5, 20);             
delta_grid = logspace(-2, -0.5, 20);              
alpha_bc   = 2e-11;                            

grids = {alpha_grid, gamma_grid, delta_grid};
keepJ = 1:opt.every:numel(tau);
pts = zeros(0, 3);
for k = 1:3
    [jj, ii] = meshgrid(keepJ, 1:opt.every:numel(grids{k}));
    pts = [pts; k*ones(numel(ii), 1), ii(:), jj(:)]; 
end
mine = pts(opt.chunk(1):opt.chunk(2):end, :);

if ~isfolder(opt.runDir), mkdir(opt.runDir); end
if license('test', 'Distrib_Computing_Toolbox') && isempty(gcp('nocreate')), parpool('Processes'); end
for q = 1:size(mine, 1)
    k = mine(q, 1);  i = mine(q, 2);  j = mine(q, 3);
    f = point_file(opt.runDir, k, i, j, opt.nReps, opt.Tend);
    if isfile(f), continue; end
    Pq = P;  Pq.alpha = alpha_bc;  Pq.m = 1/tau(j);
    switch k
        case 1, Pq.alpha = alpha_grid(i);
        case 2, Pq.gamma = gamma_grid(i);
        case 3, Pq.delta = delta_grid(i);
    end
    t0 = tic;
    elim = false(opt.nReps, 1);
    Tend = opt.Tend;
    parfor r = 1:opt.nReps
        rng(r + 1000*j + 1e5*i + 1e7*k, 'twister');
        out = simulate_niche_hybrid(Pq, Tend, dt, Lambda);
        elim(r) = out.eliminated;
    end
    nElim = nnz(elim);
    save(f, 'nElim');
    fprintf('panel %d, row %d, column %d: %d of %d eliminated (%.1f min)\n', k, i, j, nElim, opt.nReps, toc(t0)/60);
end

P_elim = cell(1, 3);
missing = 0;
for k = 1:3
    rows = 1:opt.every:numel(grids{k});
    P_elim{k} = nan(numel(rows), numel(keepJ));
    for a = 1:numel(rows)
        for b = 1:numel(keepJ)
            f = point_file(opt.runDir, k, rows(a), keepJ(b), opt.nReps, opt.Tend);
            if isfile(f)
                s = load(f);  P_elim{k}(a, b) = s.nElim/opt.nReps;
            else
                missing = missing + 1;
            end
        end
    end
end
if missing > 0
    fprintf('%d of %d grid points still to run.\n', missing, size(pts, 1));
    return;
end
P_elim_a = P_elim{1};  P_elim_b = P_elim{2};  P_elim_c = P_elim{3};
tau = tau(keepJ);
alpha_grid = alpha_grid(1:opt.every:end);  gamma_grid = gamma_grid(1:opt.every:end);
delta_grid = delta_grid(1:opt.every:end);
Tend = opt.Tend;  nReps = opt.nReps;  params = P;
method = sprintf(['Reaction-based hybrid algorithm (Supplementary Section A, dt = %g, Lambda = %g) with the ' ...
                  'niche reactions, %d realisations per grid point, seed r + 1000*column + 1e5*row + 1e7*panel ' ...
                  '(full-grid indices); elimination means B = H = 0 at day %g.'], dt, Lambda, nReps, Tend);
save(opt.outFile, 'P_elim_a', 'P_elim_b', 'P_elim_c', 'tau', 'alpha_grid', 'gamma_grid', 'delta_grid', ...
     'alpha_bc', 'Tend', 'nReps', 'params', 'method');
fprintf('Wrote %s\n', opt.outFile);
end

function f = point_file(runDir, k, i, j, nReps, Tend)
f = fullfile(runDir, sprintf('panel%d_row%02d_col%02d_n%d_T%g.mat', k, i, j, nReps, Tend));
end

function P = default_params_table1()
P.k1     = 0.2;
P.k2     = 4e-10;
P.k3     = 1e-9;
P.k4     = 0.1;
P.gamma  = 0.3;
P.delta  = 0.2;
P.eps    = 0.01;
P.K      = 1e12;
P.N      = 6;
P.Bhalf  = 1e9;

P.CAR0   = 4.1e8;
P.M0     = .468 * P.CAR0;
P.A0     = 0;
P.E0each = (.532 * P.CAR0)/P.N;

P.B0     = 2e11;   
end

function out = simulate_niche_hybrid(P, Tend, dt, Lambda)
N = P.N;
tgrid = (0:dt:Tend).';

idxB = 1;
idxE = 2:(N+1);
idxA = N + 2;
idxM = N + 3;
idxH = N + 4;
nState = N + 4;

S = beam_stoichiometry(P);

x = zeros(nState,1);
x(idxB) = P.B0;
x(idxE) = P.E0each;
x(idxA) = P.A0;
x(idxM) = P.M0;


wasFullyStochastic = false(nState,1);


odeOpts = odeset('RelTol',1e-7,'AbsTol',1e-10, ...
                  'NonNegative',1:nState);

out.eliminated = false;
out.escaped = false;
for k = 2:numel(tgrid)
    [x,wasFullyStochastic] = hybrid_reaction_interval( ...
        P, x, tgrid(k-1), tgrid(k), dt, Lambda, S, odeOpts, wasFullyStochastic);

    x = max(x,0);

    if x(idxB) == 0 && x(idxH) == 0
        out.eliminated = true;
        break;
    end
    if x(idxB) > 0.9*P.K || sum(x(idxE)) + x(idxA) + x(idxM) == 0
        out.escaped = true;
        break;
    end
end
out.x = x;
end

function [x,wasFullyStochastic] = hybrid_reaction_interval( ...
    P, x, tStart, tEnd, dtClass, Lambda, S, odeOpts, wasFullyStochastic)

nState = numel(x);
t = tStart;
timeTol = 100*eps(max(1,abs(tEnd)));

while t < tEnd - timeTol

    alpha = beam_propensities(x,P);
    detMask = classify_reactions(x,alpha,P,dtClass,Lambda);
    stochMask = ~detMask;

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

    if ~any(stochMask)
        [~,xx] = ode45(@(tt,xx) deterministic_rhs(tt,xx,P,S,detMask), ...
                       [t tEnd], x, odeOpts);
        x = xx(end,:).';
        t = tEnd;
        continue;
    end

    u = rand;
    hazardTarget = -log(max(1-u,realmin));

    z0 = [x; 0];
    eventOpts = odeset(odeOpts,'Events', ...
        @(tt,zz) hazard_event(tt,zz,hazardTarget));

    [~,zz,te,ze] = ode45( ...
        @(tt,zz) augmented_hybrid_rhs(tt,zz,P,S,detMask,stochMask,nState), ...
        [t tEnd], z0, eventOpts);

    if isempty(te)
        x = zz(end,1:nState).';
        t = tEnd;
        continue;
    end

    xEvent = ze(end,1:nState).';
    alphaEvent = beam_propensities(xEvent,P);
    rates = feasible_stochastic_rates(xEvent,alphaEvent,S,stochMask);
    totalRate = sum(rates);

    if totalRate <= 0
        x = xEvent;
        t = te(end);
        continue;
    end

    r = rand * totalRate;
    j = find(cumsum(rates) >= r,1,'first');

    x = xEvent + S(:,j);

    x(x < 0 & x > -1e-9) = 0;
    if any(x < -1e-9)
        error('Hybrid update produced a negative population after reaction %d.',j);
    end

    t = te(end);
end
end

function dx = deterministic_rhs(~,x,P,S,detMask)
alpha = beam_propensities(x,P);
if any(detMask)
    dx = S(:,detMask) * alpha(detMask);
else
    dx = zeros(size(x));
end
end


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


function [value,isterminal,direction] = hazard_event(~,z,hazardTarget)
value = z(end) - hazardTarget;
isterminal = 1;
direction = 1;
end



function alpha = beam_propensities(x,P)
N = P.N;
x = max(x,0);

B = x(1);
E = x(2:N+1);
A = x(N+2);
M = x(N+3);
H = x(N+4);
Etot = sum(E);

if P.Bhalf + B > 0
    fracEff = B/(P.Bhalf + B);
    fracMem = P.Bhalf/(P.Bhalf + B);
else
    fracEff = 0;
    fracMem = 1;
end

alpha = zeros(N+9,1);

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
alpha(N+8) = P.alpha * B;   % niche entry
alpha(N+9) = P.m * H;       % niche exit

alpha(~isfinite(alpha) | alpha < 0) = 0;
end


function S = beam_stoichiometry(P)
N = P.N;
nState = N + 4;
nRxn = N + 9;

S = zeros(nState,nRxn);

idxB = 1;
idxA = N + 2;
idxM = N + 3;

S(idxB,1) = 1;

S(idxB,2) = -1;

S(idxB,3) = -1;

S(idxM,4) = -1;
S(idxA,4) = 1;


S(idxA,5) = -1;
S(idxM,5) = 2;


S(idxA,6) = -1;
S(2,6) = 2;

for i = 1:N-1
    j = 6 + i;
    idxEi = 1 + i;
    idxEip1 = 2 + i;
    S(idxEi,j) = -1;
    S(idxEip1,j) = 2;
end

S(N+1,N+6) = -1;

S(idxM,N+7) = -1;

S(idxB,N+8) = -1;
S(N+4,N+8) = 1;

S(N+4,N+9) = -1;
S(idxB,N+9) = 1;
end


function detMask = classify_reactions(x,alpha,P,dtClass,Lambda)
N = P.N;

B = max(x(1),0);
E = max(x(2:N+1),0);
A = max(x(N+2),0);
M = max(x(N+3),0);
H = max(x(N+4),0);
Etot = sum(E);


fast = false(N+9,1);
pos = alpha > 0;
fast(pos) = (1./alpha(pos)) < dtClass;

reactantsLarge = false(N+9,1);
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
reactantsLarge(N+8) = B > Lambda;                 % niche entry
reactantsLarge(N+9) = H > Lambda;                 % niche exit

detMask = fast & reactantsLarge;
end

function n = stochastic_round_nonnegative(x)
if x <= 0
    n = 0;
    return;
end

n0 = floor(x);
p = x - n0;
n = n0 + (rand < p);
end
