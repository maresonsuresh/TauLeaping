function results = jonesadaptivetau(X, V, c, tfinal, epsilon, delta, Ncrit, Nstiff, SSAfactor, SSAsteps)
% Inputs:
% X - initial state []
% Q - stoichiometric matrix
% c - rate constants []
% tfinal - final simulation time
% Other parameters

t = 0;
type = 0;
isNegative = false;
results.time = [];
results.states = [];
results.time(1) = t;
results.states(:,1) = X;
rng(2);
tic;

while t < tfinal
    propensities = calculatePropensities(X, c);
    a0 = sum(propensities);

    if a0 == 0
        break; % No more reactions possible
    end

    % Determine step type and size
    criticalReactions = listOfCriticalReactions(X, V, propensities, Ncrit);
    [tau, type, crit] = computeTimeStep(X, V, c, propensities, criticalReactions, epsilon, delta, Nstiff);

    fprintf('t = %.4e | tau = %.4e | Type: %s | X = [%d, %d, %d, %d, %d]\n', ...
        t, tau, ...
        string(type), X(1), X(2), X(3), X(4), X(5));

    % Decide between SSA or tau-leaping
    if tau <= SSAfactor * (1.0/a0)
        % Do SSA steps
        [X, t, type] = executeSSA(X, V, c, propensities, t, tfinal, type, SSAsteps);
    else
        % Do tau-leaping
        attempts = 0;
        while attempts < 10 % Limit rejection
            if type == 0
                fire = sampling(X, V, c, propensities, tau, type, criticalReactions, crit);
            else
                fire = implicit_sampling(X, V, c, propensities, tau, criticalReactions);
            end

            % Apply reactions
            X_new = X;
            for j = 1:size(V,2)
                X_new = X_new + V(:,j) * fire(j);
            end

            % Check for negative populations
            if all(X_new >= 0)
                X = X_new;
                t = t + tau;
                isNegative = false;
                break;
            else
                tau = tau * 0.5;
                isNegative = true;
                attempts = attempts + 1;
            end
        end
    end

    % Store results
    results.time(end+1) = t;
    results.states(:,end+1) = X;
end

elapsed_time = toc;
fprintf('Simulation time: %.4f seconds\n', elapsed_time);
figure;
hold on;

% PLOTTING!!!!
%plot(results.time, results.states(1,:), 'b', 'LineWidth', 2, 'DisplayName', 'A');
plot(results.time, results.states(2,:), 'r', 'LineWidth', 2, 'DisplayName', 'CC');
plot(results.time, results.states(3,:), 'g', 'LineWidth', 2, 'DisplayName', 'T');
plot(results.time, results.states(4,:), 'm', 'LineWidth', 2, 'DisplayName', 'T*');
plot(results.time, results.states(5,:), 'k', 'LineWidth', 2, 'DisplayName', 'V');
xlabel('Time (days)');
ylabel('Molecule Count');
legend('Location', 'best');
grid on;
xlim([0, tfinal]);
end


function critical_reactions = listOfCriticalReactions(X, V, propensities, Ncrit)

num_reactions = size(V, 2);
critical_reactions = zeros(1, num_reactions); % Preallocate
crit_count = 0; % Counter for crit

for ir = 1:num_reactions
    if propensities(ir) > 0
        lj = inf;
        for is = 1:size(V, 1)
            nu = V(is, ir); % Stoich
            if nu < 0
                lj_candidate = floor(X(is) / abs(nu));
                if lj_candidate < lj
                    lj = lj_candidate;
                end
            end
        end
        % If drop below Ncrit, critical
        if lj < Ncrit
            crit_count = crit_count + 1;
            critical_reactions(crit_count) = ir;
        end
    end
end
% Trim the preallocate
critical_reactions = critical_reactions(1:crit_count);
end

function [tau, type, crit] = computeTimeStep(X, V, ~, propensities, criticalReactions, epsilon, delta, Nstiff)
tauPrimeExp = inf;  
tauPrimeImp = inf; 
tau1 = inf;        
tau2 = inf;        
rev = [10, 8, 2, 4]; % Indices of reversible reaction pairs (1-based indexing)
non_critical = [];
in_pec = [];  
num_species = size(V, 1);
num_reactions = size(V, 2);

muHat = zeros(num_species, 1);
sigmaHat2 = zeros(num_species, 1);
varHat = zeros(num_species, 1);

% Compute hor and nuHor
% Highest order of reaction for each species
hor = [2;  
       1;  
       2;  
       1; 
       2];
nuHor = [1;  
         1; 
         1;
         1;  
         1];
% 1. STEP: EXPLICIT TAU
if isempty(criticalReactions)
    non_critical = 1:num_reactions;
else
    non_critical = setdiff(1:num_reactions, criticalReactions);
end

% Compute muHat and sigmaHat2 for non-critical reactions
[muHat, sigmaHat2] = computeMuHatSigmaHat2(X, V, propensities, non_critical);

% Compute explicit tau
tau = inf;
a0 = sum(propensities);

for is = 1:num_species
    varHat(is) = sigmaHat2(is) - (1.0/a0) * muHat(is) * muHat(is);
end

for is = 1:num_species
    xi = X(is);
    switch hor(is)
        case 0
            continue;
        case 1
            epsi = epsilon;
            epsixi = epsi * xi;
            epsixi = max(epsixi, 1.0);
            tau = min(tau, epsixi/abs(muHat(is)));
            epsixisq = epsixi^2;
            tau = min(tau, epsixisq/varHat(is));
        case 2
            if nuHor(is) == 1
                epsi = 0.5*epsilon;
            else
                epsi = epsilon*(xi-1.0)/(2.0*(xi-1.0)+1.0);
            end
            epsixi = epsi * xi;
            epsixi = max(epsixi, 1.0);
            tau = min(tau, epsixi/abs(muHat(is)));
            epsixisq = epsixi^2;
            tau = min(tau, epsixisq/varHat(is));
        case 3
            if nuHor(is) == 1
                epsi = 0.3333333333*epsilon;
            elseif nuHor(is) == 2
                epsi = epsilon*(xi-1)/(3.0*(xi-1)+1.5);
            else
                epsi = epsilon*(xi-1)*(xi-2)/(3.0*(xi-1)*(xi-2)+(xi-2)+2.0*(xi-1));
            end
            epsixi = epsi * xi;
            epsixi = max(epsixi, 1.0);
            tau = min(tau, epsixi/abs(muHat(is)));
            epsixisq = epsixi^2;
            tau = min(tau, epsixisq/varHat(is));
    end
end
tauPrimeExp = tau;

% 2. STEP: IMPLICIT TAU
if isempty(rev)
    tauPrimeImp = tauPrimeExp;
else
    % Check for reactions in partial equilibrium
    for i = 1:2:length(rev)
        temp1 = propensities(rev(i));   
        temp2 = propensities(rev(i+1)); 

        if abs(temp1 - temp2) < delta * abs(temp1 + temp2)
            in_pec = [in_pec, rev(i), rev(i+1)];
        end
    end

    if isempty(in_pec)
        tauPrimeImp = tauPrimeExp;
    else
        % Update non_critical by removing PEC
        non_critical = setdiff(non_critical, in_pec);
        % Recompute muHat and sigmaHat2
        [muHat, sigmaHat2] = computeMuHatSigmaHat2(X, V, propensities, non_critical);
        % Recompute tau with updated muHat and sigmaHat2
        tau = inf;
        a0 = sum(propensities);

        for is = 1:num_species
            varHat(is) = sigmaHat2(is) - (1.0/a0) * muHat(is) * muHat(is);
        end

        for is = 1:num_species
            xi = X(is);
            switch hor(is)
                case 0
                    continue;
                case 1
                    epsi = epsilon;
                    epsixi = epsi * xi;
                    epsixi = max(epsixi, 1.0);
                    tau = min(tau, epsixi/abs(muHat(is)));
                    epsixisq = epsixi^2;
                    tau = min(tau, epsixisq/varHat(is));
                case 2
                    if nuHor(is) == 1
                        epsi = 0.5*epsilon;
                    else
                        epsi = epsilon*(xi-1.0)/(2.0*(xi-1.0)+1.0);
                    end
                    epsixi = epsi * xi;
                    epsixi = max(epsixi, 1.0);
                    tau = min(tau, epsixi/abs(muHat(is)));
                    epsixisq = epsixi^2;
                    tau = min(tau, epsixisq/varHat(is));
                case 3
                    if nuHor(is) == 1
                        epsi = 0.3333333333*epsilon;
                    elseif nuHor(is) == 2
                        epsi = epsilon*(xi-1)/(3.0*(xi-1)+1.5);
                    else
                        epsi = epsilon*(xi-1)*(xi-2)/(3.0*(xi-1)*(xi-2)+(xi-2)+2.0*(xi-1));
                    end
                    epsixi = epsi * xi;
                    epsixi = max(epsixi, 1.0);
                    tau = min(tau, epsixi/abs(muHat(is)));
                    epsixisq = epsixi^2;
                    tau = min(tau, epsixisq/varHat(is));
            end
        end
        tauPrimeImp = tau;
    end
end

% 3. Determine tau
if tauPrimeImp > Nstiff * tauPrimeExp
    type = 1;       % Stiff system - use implicit tau
    tau1 = tauPrimeImp;
else
    type = 0;       % Non-stiff system - use explicit tau
    tau1 = tauPrimeExp;
end
% 4. Compute time of next critical reaction
if isempty(criticalReactions)
    tau2 = inf;
else
    ac0 = sum(propensities(criticalReactions));
    if ac0 > 0
        tau2 = exprnd(1/ac0);
    else
        tau2 = inf;
    end
end
if tau1 < tau2
    tau = tau1;
    crit = 0;
else
    tau = tau2;
    crit = 1;

    if (type == 0) || (type == 1 && tau2 < tauPrimeExp)
        type = 0;
    else
        type = 1;
    end
end
end

function [muHat, sigmaHat2] = computeMuHatSigmaHat2(~, V, propensities, non_critical)
num_species = size(V, 1);
muHat = zeros(num_species, 1);
sigmaHat2 = zeros(num_species, 1);

% Consider only non-critical reactions
for i = 1:length(non_critical)
    ir = non_critical(i);
    ri_propensity = propensities(ir);
    nuChanges = V(:, ir);

    for is = 1:num_species
        tmpfloat = nuChanges(is) * ri_propensity;
        muHat(is) = muHat(is) + tmpfloat;
        sigmaHat2(is) = sigmaHat2(is) + nuChanges(is) * tmpfloat;
    end
end
end

function [X, t, type] = executeSSA(X, V, c, ~, t, tfinal, type, SSAsteps)
% Determine number of SSA steps based on previous step type
if type == 0
    steps = SSAsteps; 
    fprintf('Executing %d SSA steps\n', SSAsteps);
else
    steps = SSAsteps/10; 
    fprintf('Executing %d SSA steps\n', SSAsteps/10);
end
type = 0; % Always return type 0

count = 0;
while count < steps && t < tfinal
    % Compute propensities
    propensities = calculatePropensities(X, c);
    a0 = sum(propensities);

    if a0 == 0
        t = inf;
        break;
    end

    % Generate time step
    dt = exprnd(1/a0);
    
    % Select reaction
    r1 = rand() * a0;
    cumulative = 0;
    reactionIndex = 0;
    
    for j = 1:length(propensities)
        cumulative = cumulative + propensities(j);
        if cumulative > r1
            reactionIndex = j;
            break;
        end
    end

    if reactionIndex > 0
        X = X + V(:, reactionIndex);
        X = max(X, 0);  % Prevent negative
        t = t + dt;
    else
        t = inf;
        break;
    end
    count = count + 1;
end
end

function fire = sampling(~, ~, ~, propensities, tau, type, criticalReactions, crit)
% Initialize fire vector for 11 reactions
fire = zeros(1, 11);

% Debug: Verify inputs
if length(propensities) ~= 11
    error('Propensities must have 11 elements, got %d', length(propensities));
end
if ~isempty(criticalReactions) && any(criticalReactions < 1 | criticalReactions > 11)
    error('Critical reactions must be between 1-11, got %s', mat2str(criticalReactions));
end

% Explicit sampling for non-critical reactions
if type == 0
    for j = 1:11 % Loop through all 11 reactions
        if isempty(criticalReactions) || ~ismember(j, criticalReactions)
            aj = propensities(j);
            if aj > 0
                fire(j) = poissrnd(aj * tau);
            end
        end
    end
else
    error('Implicit sampling should use implicit_sampling() function');
end

% Handle critical reaction firing
if crit == 1 && ~isempty(criticalReactions)
    ac0 = sum(propensities(criticalReactions));
    if ac0 > 0
        r1 = rand() * ac0;
        cumulative = 0;
        for k = 1:length(criticalReactions)
            idx = criticalReactions(k);
            if idx < 1 || idx > 11
                error('Invalid reaction index %d in criticalReactions', idx);
            end
            cumulative = cumulative + propensities(idx);
            if cumulative > r1
                fire(idx) = fire(idx) + 1;
                break;
            end
        end
    end
end

% Final verification
if length(fire) ~= 11
    error('Fire vector has %d elements (expected 11)', length(fire));
end
end

%NEW VERSION USING NEWTON RAPHSON
function fire = implicit_sampling(X, V, c, propensities, tau, critical)
% implicit tau-leaping with Newton-Raphson for Jones model
% Compute Poisson variates
k = poissrnd(propensities * tau);
a_tau = propensities * tau;
B = X + V * (k - a_tau);

% Newton-Raphson solver
max_iter = 100;
tol = 1e-8;
damp_factor = 0.5;
min_step_size = 1e-10;

% Initialize variables
roots = X;              % Initial guess
converged = false;
iter = 0;

% Newton-Raphson iteration
while ~converged && iter < max_iter
    iter = iter + 1;
    % Calculate propensities and Jacobian
    [a_current, J] = calculatePropensitiesWithJacobian(roots, c, V, tau);
    % Compute function value
    F = roots - B - V*a_current*tau;
    % Check convergence
    if norm(F) < tol
        converged = true;
        break;
    end
    
    % Compute Newton step: J*dx = -F
    dx = J \ (-F);
    step_size = damp_factor;

    % Ensure non-negative populations
    while step_size >= min_step_size
        new_roots = roots + step_size * dx;
        if all(new_roots >= 0)
            roots = new_roots;
            break;
        else
            step_size = step_size * damp_factor;
        end
    end

    % Fallback to explicit approximation if needed
    if step_size < min_step_size
        warning('Newton-Raphson failed - using explicit approximation');
        roots = X + V*(propensities*tau);
        break;
    end
end

if ~converged && iter == max_iter
    warning('Newton-Raphson reached max iterations');
end

% Compute reaction counts
a_roots = calculatePropensities(roots, c);
fire = round(a_roots * tau + k - a_tau);
fire = max(fire, 0);  % Ensure non-negativity

% UPDATED DEBUG PRINT FOR 5 SPECIES
fprintf('[Implicit] X_old = [%d %d %d %d %d] -> X_new = [%d %d %d %d %d] (NR iters: %d)\n', ...
    round(X(1)), round(X(2)), round(X(3)), round(X(4)), round(X(5)), ...
    round(roots(1)), round(roots(2)), round(roots(3)), round(roots(4)), round(roots(5)), iter);

% Zero critical reactions
fire(critical) = 0;
end

function [a, J] = calculatePropensitiesWithJacobian(X, c, V, tau)
% Jones model version
num_species = length(X); % 5 species
num_reactions = 11;      % 11 reactions
a = zeros(num_reactions, 1);

% Calculate propensities (same as calculatePropensities)
a = [
    c(1)*X(1)*X(3);                     % v1
    c(7)*c(11)*X(4);                     % v10
    c(8)*c(12)*X(2);                     % v11
    c(2);                                % v2
    c(3)*X(1)*X(3)/(c(4) + X(1));        % v3
    c(5)*X(3)*X(5);                      % v4
    c(6)*X(3);                           % v5
    c(7)*X(4);                           % v6
    c(8)*X(2);                           % v7
    c(10)*X(5);                          % v8
    c(9)*c(5)*X(3)*X(5)                  % v9
    ];

% Jacobian components (da_i/dX_j)
J_a = zeros(num_reactions, num_species);

% Row 1: v1 = c1*A*T
J_a(1,1) = c(1)*X(3);        % d(v1)/dA
J_a(1,3) = c(1)*X(1);        % d(v1)/dT
% Row 2: v10 = c7*c11*T*
J_a(2,4) = c(7)*c(11);       % d(v10)/dT*
% Row 3: v11 = c8*c12*CC
J_a(3,2) = c(8)*c(12);       % d(v11)/dCC
% Row 5: v3 = c3*A*T/(c4 + A)
J_a(5,1) = c(3)*X(3)*c(4)/(c(4) + X(1))^2; % d(v3)/dA
J_a(5,3) = c(3)*X(1)/(c(4) + X(1));        % d(v3)/dT
% Row 6: v4 = c5*T*V
J_a(6,3) = c(5)*X(5);        % d(v4)/dT
J_a(6,5) = c(5)*X(3);        % d(v4)/dV
% Row 7: v5 = c6*T
J_a(7,3) = c(6);             % d(v5)/dT
% Row 8: v6 = c7*T*
J_a(8,4) = c(7);             % d(v6)/dT*
% Row 9: v7 = c8*CC
J_a(9,2) = c(8);             % d(v7)/dCC
% Row 10: v8 = c10*V
J_a(10,5) = c(10);           % d(v8)/dV
% Row 11: v9 = c9*c5*T*V
J_a(11,3) = c(9)*c(5)*X(5);  % d(v9)/dT
J_a(11,5) = c(9)*c(5)*X(3);  % d(v9)/dV

% Compute full Jacobian: J = I - tau*V*J_a
J = eye(num_species) - tau * V * J_a;
end

function [a, a0] = calculatePropensities(X, c)
% X = [A; CC; T; T*; V]
a = [
    c(1)*X(1)*X(3);                    
    c(7)*c(11)*X(4);                     
    c(8)*c(12)*X(2);                    
    c(2);                               
    c(3)*X(1)*X(3)/(c(4) + X(1));       
    c(5)*X(3)*X(5);                     
    c(6)*X(3);                        
    c(7)*X(4);                        
    c(8)*X(2);                         
    c(10)*X(5);                  
    c(9)*c(5)*X(3)*X(5)         
    ];
a = max(a, 0);  % Ensure non-negative propensities
a0 = sum(a);    % Total propensity
end