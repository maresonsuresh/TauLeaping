function DDM_ATL()
    % Parameters
    c = [0.05, 50, 1e6, 0.05];      % Rate constants
    V = [-1  -2   2    0;           % Stoichiometric matrix (S1, S2, S3)
          0   1  -1   -1;
          0   0   0    1];
    X = [400; 800; 0];              % Initial state
    t = 0; tfinal = 1;              % Time parameters
    Nc = 10;                        % Critical reaction threshold
    epsilon = 0.03;                 % Error control parameter
    max_steps = 1e5;                % Maximum simulation steps
    results = zeros(max_steps,4);   % Storage for results
    results(1,:) = [t X'];          % Initial state
    ii = 1; rng('shuffle'); tic;    % Initialize
    
    % Main simulation loop
    while t < tfinal
        % Step 3: Compute propensities
        a = [c(1)*X(1);
             c(2)*X(1)*(X(1)-1)/2;
             c(3)*X(2);
             c(4)*X(2)];
        a0 = sum(a);
        if a0 == 0, break; end

        % Step 4: Identify critical reactions
        critical = false(4,1);
        for j = 1:4
            if a(j) > 0
                reactants = find(V(:,j) < 0);
                L_j = min(floor(X(reactants)./abs(V(reactants,j))));
                critical(j) = (L_j <= Nc);
            end
        end

        % Step 5: Compute tau candidates
        [tau_ex, tau_im] = compute_taus(X, a, critical, V, epsilon);
        
        % Steps 6-10: Determine stiffness
        if tau_im > 100*tau_ex
            stiff = true;
            tau1 = tau_im;
        else
            stiff = false;
            tau1 = tau_ex;
        end

        % Steps 11-12: Check if SSA should be used
        if tau1 <= 10/a0
            [t, X, ii] = ssa_steps(t, X, a, V, ii, results, tfinal);
            continue;
        end

        % Steps 14-32: Generate firings
        crit_idx = find(critical);
        if ~isempty(crit_idx)
            a0_c = sum(a(crit_idx));
            tau2 = exprnd(1/a0_c);
        else
            tau2 = inf;
        end
        
        k = zeros(4,1);
        if tau2 > tau1
            tau = tau1;
            if stiff
                k(~critical) = implicit_firings(X, a, V, find(~critical), tau);
            else
                k(~critical) = poissrnd(a(~critical)*tau);
            end
        else
            tau = tau2;
            jc = crit_idx(find(cumsum(a(crit_idx)) >= rand*a0_c, 1));
            k(jc) = 1;
            if tau2 < tau_ex || ~stiff
                k(~critical) = poissrnd(a(~critical)*tau);
            else
                k(~critical) = implicit_firings(X, a, V, find(~critical), tau);
            end
        end

        % Steps 33-39: Validate and update
        X_new = X + V*k;
        if any(X_new < 0)
            tau1 = tau1/2;
            continue;
        else
            X = X_new;
            t = t + tau;
            ii = ii + 1;
            results(ii,:) = [t X'];
        end
    end

    % Post-processing
    results = results(1:ii,:);
    fprintf('Simulation completed in %.2f seconds\n', toc);
    
    % Plotting
    figure('Color','white');
    plot(results(:,1), results(:,2), 'LineWidth',2, 'DisplayName','S1');
    hold on;
    plot(results(:,1), results(:,3), 'LineWidth',2, 'DisplayName','S2');
    plot(results(:,1), results(:,4), 'LineWidth',2, 'DisplayName','S3');
    hold off;
    xlabel('Time (s)'); ylabel('Molecule Count');
    legend('Location','best'); grid on;
    set(gca, 'FontSize',12);
end

%% Helper functions
function [tau_ex, tau_im] = compute_taus(X, a, critical, V, eps)
    % Explicit tau (non-critical reactions)
    Jncr = find(~critical);
    if isempty(Jncr)
        mu_ex = 0; sigma2_ex = 0;
    else
        mu_ex = sum(V(:,Jncr).*a(Jncr)', 2);
        sigma2_ex = sum(V(:,Jncr).^2.*a(Jncr)', 2);
    end
    
    % Implicit tau (non-critical reactions)
    % Note: Original pseudocode doesn't include equilibrium checks
    mu_im = mu_ex;
    sigma2_im = sigma2_ex;
    
    % Tau calculation
    g_i = ones(3,1);
    tau_ex = min(max(eps*X./g_i,1)./abs(mu_ex), max(eps*X./g_i,1).^2./sigma2_ex);
    tau_im = min(max(eps*X./g_i,1)./abs(mu_im), max(eps*X./g_i,1).^2./sigma2_im);
    tau_ex = min(tau_ex); tau_im = min(tau_im);
end

function k = implicit_firings(~, a, ~, reactions, tau)
    % Simplified implicit firing calculation
    % Note: In practice, you'd use a more sophisticated implicit method
    k = zeros(4,1);
    if ~isempty(reactions)
        % Simple approximation - in real implementation use proper implicit solver
        k(reactions) = poissrnd(a(reactions).*tau);
    end
end

function [t_new, X_new, ii_new] = ssa_steps(t, X, a, V, ii, results, tfinal)
    for ssa = 1:100
        dt = log(1/rand)/sum(a);
        j = find(cumsum(a) >= rand*sum(a),1);
        X = X + V(:,j);
        t = t + dt;
        ii = ii + 1;
        results(ii,:) = [t X'];
        if t >= tfinal, break; end
    end
    t_new = t; X_new = X; ii_new = ii;
end