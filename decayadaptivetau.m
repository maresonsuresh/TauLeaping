function results = decayadaptivetau(X, V, c, tfinal, epsilon, delta, Ncrit, Nstiff, SSAfactor, SSAsteps)
% Inputs:
% X - initial state [S1; S2; S3]
% V - stoichiometric matrix
% c - rate constants [C1, C2, C3, C4]
% tfinal - final simulation time
% Other parameters as defined in your original input

% Initialize
t = 0;
type = 0;
isNegative = false;
results.time = [];
results.states = [];

% Store initial condition
results.time(1) = t;
results.states(:,1) = X;
rng(2);
tic;
while t < tfinal
    % Compute propensities
    propensities = calculatePropensities(X, c);
    a0 = sum(propensities);

    if a0 == 0
        break; % No more reactions possible
    end

    % Determine step type and size
    criticalReactions = listOfCriticalReactions(X, V, propensities, Ncrit);
    [tau, type, crit] = computeTimeStep(X, V, c, propensities, criticalReactions, epsilon, delta, Nstiff);

    fprintf('t = %.4e | tau = %.4e (size: %.4e) | Type: %s | X = [%d, %d, %d]\n', ...
    t, tau, tau, ...
    string(type), X(1), X(2), X(3)); 

    % Decide between SSA or tau-leaping
    if tau <= SSAfactor * (1.0/a0)
        % Do SSA steps
        [X, t, type] = executeSSA(X, V, c, propensities, t, tfinal, type, SSAsteps);
    else
        % Do tau-leaping
        attempts = 0;
        while attempts < 10 % Limit rejection attempts
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
% Plot all species on one graph
figure;
hold on;

% Plot each species with different colors and styles
plot(results.time, results.states(1,:), 'b', 'LineWidth', 2, 'DisplayName', 'S1');
plot(results.time, results.states(2,:), 'r', 'LineWidth', 2, 'DisplayName', 'S2'); 
plot(results.time, results.states(3,:), 'g', 'LineWidth', 2, 'DisplayName', 'S3');


% Add labels and legend
xlabel('Time (seconds)');
ylabel('Molecule Count');
legend('Location', 'best');
grid on;

% Adjust axes for better visibility
xlim([0, tfinal]);
ylim([0 max(results.states(:))*1.1]); % 10% padding
end


function critical_reactions = listOfCriticalReactions(X, V, propensities, Ncrit)
% X: current state vector [S1; S2; S3]
% V: stoichiometric matrix (species x reactions)
% propensities: vector of current reaction propensities
% Ncrit: critical threshold for species count

num_reactions = size(V, 2);
critical_reactions = zeros(1, num_reactions); % Preallocate maximum possible size
crit_count = 0; % Counter for actual critical reactions

for ir = 1:num_reactions
    if propensities(ir) > 0
        lj = inf; % Initialize with infinity

        % Check each species affected by this reaction
        for is = 1:size(V, 1)
            nu = V(is, ir); % Stoichiometric coefficient

            % Only consider negative changes (consumption)
            if nu < 0
                lj_candidate = floor(X(is) / abs(nu));
                if lj_candidate < lj
                    lj = lj_candidate;
                end
            end
        end

        % If any species would drop below Ncrit, mark as critical
        if lj < Ncrit
            crit_count = crit_count + 1;
            critical_reactions(crit_count) = ir;
        end
    end
end
% Trim the preallocated array to actual size
critical_reactions = critical_reactions(1:crit_count);
end

function [tau, type, crit] = computeTimeStep(X, V, ~, propensities, criticalReactions, epsilon, delta, Nstiff)
% Inputs:
% X - current state vector [S1; S2; S3]
% V - stoichiometric matrix (species x reactions)
% c - reaction rate constants [C1, C2, C3, C4]
% propensities - current reaction propensities
% criticalReactions - list of critical reaction indices
% epsilon, delta, Nstiff - algorithm parameters

% Outputs:
% tau - computed time step
% type - 0 for explicit, 1 for implicit
% crit - number of critical reactions to fire

% 1. Initialize variables
tauPrimeExp = inf;  % explicit tau
tauPrimeImp = inf;  % implicit tau
tau1 = inf;         % tau for non-critical reactions
tau2 = inf;         % time of 1st critical reaction

% For your specific dimerization model, reversible reactions are R2 and R3
rev = [2, 3]; % Indices of reversible reaction pairs (1-based indexing)

non_critical = []; % list of reactions used to compute muHat, sigmaHat2
in_pec = [];       % list of reactions in partial equilibrium

num_species = size(V, 1);
num_reactions = size(V, 2);

% Initialize arrays
muHat = zeros(num_species, 1);
sigmaHat2 = zeros(num_species, 1);
varHat = zeros(num_species, 1);

% Compute hor and nuHor (placeholder - we'll implement this properly later)
% For now, we'll assume hor=1 for all species as a simplification
hor = ones(num_species, 1);
nuHor = ones(num_species, 1);

% 1. STEP: EXPLICIT TAU
%-----------------------
if isempty(criticalReactions)
    non_critical = 1:num_reactions;
else
    % Create list of non-critical reactions
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
%-----------------------
if isempty(rev)
    tauPrimeImp = tauPrimeExp;
else
    % Check for reactions in partial equilibrium
    for i = 1:2:length(rev)
        temp1 = propensities(rev(i));    % forward reaction
        temp2 = propensities(rev(i+1));  % backward reaction

        if abs(temp1 - temp2) < delta * abs(temp1 + temp2)
            in_pec = [in_pec, rev(i), rev(i+1)];
        end
    end

    if isempty(in_pec)
        tauPrimeImp = tauPrimeExp;
    else
        % Update non_critical list by removing PEC reactions
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

% 3. Determine tau and stiffness type
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
        tau2 = exprnd(1/ac0); % Exponential random variate
    else
        tau2 = inf;
    end
end

% 5. Choose final tau and number of critical reactions
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
% Inputs:
% X - current state vector [S1; S2; S3]
% V - stoichiometric matrix (species x reactions)
% propensities - current reaction propensities
% non_critical - list of non-critical reaction indices

% Outputs:
% muHat - vector of muHat values for each species
% sigmaHat2 - vector of sigmaHat2 values for each species

num_species = size(V, 1);
muHat = zeros(num_species, 1);
sigmaHat2 = zeros(num_species, 1);

% Consider only non-critical reactions
for i = 1:length(non_critical)
    ir = non_critical(i);
    ri_propensity = propensities(ir);

    % Get changes for this reaction (columns of V)
    nuChanges = V(:, ir);

    for is = 1:num_species
        tmpfloat = nuChanges(is) * ri_propensity;
        muHat(is) = muHat(is) + tmpfloat;
        sigmaHat2(is) = sigmaHat2(is) + nuChanges(is) * tmpfloat;
    end
end
end

function [X, t, type] = executeSSA(X, V, c, ~, t, tfinal, type, SSAsteps)

% Inputs:
% X - current state vector [S1; S2; S3]
% V - stoichiometric matrix
% c - reaction rate constants [C1, C2, C3, C4]
% propensities - current reaction propensities
% t - current time
% tfinal - final simulation time
% type - 0 for explicit/SSA, 1 for implicit
% SSAsteps - number of SSA steps to execute (typically 10 or 100)

% Outputs:
% X - updated state vector
% t - updated time
% type - always returns 0 (since we're doing SSA)

% Determine number of SSA steps based on previous step type
if type == 0
    steps = SSAsteps; % Typically 100
    fprintf('SSA steps! 100\n');
else
    steps = SSAsteps/10; % Typically 10
    fprintf('SSA steps! 10\n');
end
type = 0; % Always return type 0 after SSA

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
    cummulative = 0;
    reactionIndex = 0;

    for j = 1:length(propensities)
        cummulative = cummulative + propensities(j);
        if cummulative > r1
            reactionIndex = j;
            break;
        end
    end

    if reactionIndex > 0
        % Fire the reaction
        X = X + V(:, reactionIndex);
        t = t + dt;
    else
        t = inf;
        break;
    end

    count = count + 1;
end
end

function fire = sampling(~, ~, ~, propensities, tau, type, criticalReactions, crit)
% Initialize fire vector with exactly 4 elements (one per reaction)
fire = zeros(1, 4); % 4 reactions in your system

% Debug: Verify inputs
if length(propensities) ~= 4
    error('Propensities must have 4 elements, got %d', length(propensities));
end
if any(criticalReactions < 1 | criticalReactions > 4)
    error('Critical reactions must be between 1-4, got %s', mat2str(criticalReactions));
end

% Explicit sampling
if type == 0
    for j = 1:4 % Explicitly loop through all 4 reactions
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

% Add critical reaction if needed
if crit == 1 && ~isempty(criticalReactions)
    ac0 = sum(propensities(criticalReactions));
    if ac0 > 0
        r1 = rand() * ac0;
        cumulative = 0;
        for k = 1:length(criticalReactions)
            idx = criticalReactions(k);
            if idx < 1 || idx > 4
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
if length(fire) ~= 4
    error('Programming error: fire vector has %d elements', length(fire));
end
end


%NEW VERSION USING NEWTON RAPHSON
function fire = implicit_sampling(X, V, c, propensities, tau, critical)
    % implicit tau-leaping with Newton-Raphson
    
    %compute Poisson
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
    
    %Newton-Raphson iteration
    while ~converged && iter < max_iter
        iter = iter + 1;
        
        % Calculate propensities and Jacobian
        [a_current, J] = calculatePropensitiesWithJacobian(roots, c, V, tau);
        % Computefunction value
        F = roots - B - V*a_current*tau;
        % Check convergence
        if norm(F) < tol
            converged = true;
            break;
        end
        % Compute the Newton step: J*dx = -F
        dx = J \ (-F);
        step_size = damp_factor;
        
        % Ensure populations not negative
        while step_size >= min_step_size
            new_roots = roots + step_size * dx;
            if all(new_roots >= 0)
                roots = new_roots;
                break;
            else
                step_size = step_size * damp_factor;
            end
        end
        
        % If couldn't find valid step! -> use explicit approximation
        if step_size < min_step_size
            warning('Newton-Raphson failed to converge - using explicit approximation');
            roots = X + V*(propensities*tau);
            break;
        end
    end
    
    if ~converged && iter == max_iter
        warning('Newton-Raphson reached maximum iterations - using current estimate');
    end
    
    %compute reaction counts
    a_roots = calculatePropensities(roots, c);
    fire = round(a_roots * tau + k - a_tau);
    fire = max(fire, 0);  %nonnegativity
    
    %Debug
    fprintf('[Implicit] X_old = [%d %d %d] -> X_new = [%d %d %d] (NR iters: %d)\n', ...
            round(X(1)), round(X(2)), round(X(3)), ...
            round(roots(1)), round(roots(2)), round(roots(3)), iter);
    % Zero critical reactions
    fire(critical) = 0;
end

function [a, J] = calculatePropensitiesWithJacobian(X, c, V, tau)
    % Calculate propensities and Jacobian
    % Initialize
    num_species = length(X);
    num_reactions = 4; % For this specific system
    a = zeros(num_reactions, 1);
    J = zeros(num_species, num_species); % Jacobian of F(x) = x - B - V*a(x)*tau
    
    %Calculate propensities

    a(1) = c(1)*X(1);                     
    a(2) = c(2)*X(1)*(X(1)-1)/2;     
    a(3) = c(3)*X(2);             
    a(4) = c(4)*X(2);     

    J_a = zeros(num_reactions, num_species);
    J_a(1,1) = c(1);
    J_a(2,1) = c(2)*(X(1)-0.5);
    J_a(3,2) = c(3);
    J_a(4,2) = c(4);
    %compute the full Jacobian: J = I - tau*V*J_a
    J = eye(num_species) - tau * V * J_a;
end

% OLD VERSION USING FSOLVE!
% function fire = implicit_sampling(X, V, c, propensities, tau, critical)
%     % Robust implicit tau-leaping with state update validation
% 
%     % --- Step 1: Precompute Poisson counts and B vector ---
%     k = poissrnd(propensities * tau);
%     a_tau = propensities * tau;
%     B = X + V * (k - a_tau);
% 
%     % --- Step 2: Configure fsolve with stricter settings ---
%     options = optimoptions('fsolve', ...
%         'Display', 'none', ...
%         'FunctionTolerance', 1e-8, ...
%         'StepTolerance', 1e-8, ...
%         'MaxIterations', 100000);
% 
%     % --- Step 3: Solve X = B + V*a(X)*tau ---
%     try
%         [roots, ~, exitflag] = fsolve(@(x) x - B - V*calculatePropensities(x,c)*tau, X, options);
% 
%         % Fallback if fsolve fails or returns invalid states
%         if exitflag <= 0 || any(roots < 0)
%             warning('Implicit solver failed - using explicit approximation');
%             roots = X + V*(propensities*tau);  % Explicit Euler fallback
%         end
%     catch
%         roots = X;  % Ultimate fallback
%     end
% 
%     % --- Step 4: Compute validated reaction counts ---
%     a_roots = calculatePropensities(roots, c);
%     fire = round(a_roots * tau + k - a_tau);
%     fire = max(fire, 0);  % Ensure non-negativity
% 
%     % --- Step 5: Debug output ---
%     fprintf('[Implicit] X_old = [%d %d %d] -> X_new = [%d %d %d]\n', ...
%             round(X(1)), round(X(2)), round(X(3)), ...
%             round(roots(1)), round(roots(2)), round(roots(3)));
% 
%     % Zero critical reactions
%     fire(critical) = 0;
% end

function [a, a0] = calculatePropensities(X, c)
a = [c(1)*X(1);                     % R1: S1 → ∅
    c(2)*X(1)*(X(1)-1)/2;          % R2: 2S1 → S2 (dimerization)
    c(3)*X(2);                     % R3: S2 → 2S1 (dissociation)
    c(4)*X(2)];                    % R4: S2 → S3 (decay)
a = max(a, 0);
a0 = sum(a);
end