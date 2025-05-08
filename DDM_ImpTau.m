function DDM_ImpTau()
    % Rate constants and Stoichiometric matrix
    c = [1, 10, 1000, 0.1];
    V = [-1  -2  +2   0;   % S1
          0  +1  -1  -1;   % S2
          0   0   0  +1];  % S3
    X = [400; 798; 0];     % Initial state
    t = 0;
    tfinal = 0.2;
    epsilon = 0.05;         % Tolerance for leap condition
    delta_partial = 0.05;   % Partial equilibrium threshold
    %nc = 5;                 % Critical reaction parameter
    max_newton_iter = 10;   % Max Newton iterations
    newton_tol = 1e-6;      % Newton convergence tolerance

    % Preallocate results
    max_steps = 1000;
    results = NaN(max_steps, 4);
    ii = 1;
    results(ii,:) = [t X'];
    rng(45);  % Seed for reproducibility

    tic;
    while t < tfinal
        x = X;
        a = compute_propensities(x, c);
        
        % Check partial equilibrium for R2 and R3
        if abs(a(2)-a(3)) > delta_partial * min(a(2),a(3))
            warning('R2/R3 not in partial equilibrium.');
        end
        
        % Stepsize selection using leap conditions
        [tau, critical] = select_tau(x, c, epsilon, delta_partial, tfinal-t);
        if tau <= 0
            warning('Invalid tau computed.');
            break;
        end
        
        % Generate Poisson numbers for stochastic part
        P = poissrnd(a * tau);
        
        % Compute RHS of implicit equation: x + sum(nu*(P - a*tau))
        RHS = x + V * (P - a * tau);
        
        % Newton-Raphson to solve X_new = RHS + V*(a_new*tau)
        X_guess = x; % Initial guess
        converged = false;
        for iter = 1:max_newton_iter
            a_new = compute_propensities(X_guess, c);
            G = X_guess - V * (a_new * tau) - RHS;
            if norm(G) < newton_tol
                converged = true;
                break;
            end
            % Compute Jacobian matrix
            J = compute_jacobian(X_guess, c, V, tau);
            % Update using Newton step
            delta = J \ (-G);
            X_guess = X_guess + delta;
            % Ensure non-negative populations
            X_guess = max(X_guess, 0);
        end
        
        if ~converged
            warning('Newton did not converge');
        end
        
        % Update state and time
        X = X_guess;
        t = t + tau;
        
        % Store results
        ii = ii + 1;
        if ii > max_steps
            warning('Max steps reached');
            break;
        end
        results(ii,:) = [t X'];
    end
    
    elapsed_time = toc;
    fprintf('Simulation completed in %.4f seconds with %d steps.\n', elapsed_time, ii);
    
    % Trim and plot results
    results = results(1:ii, :);
    plot_results(results);
end

function a = compute_propensities(x, c)
    a = zeros(4,1);
    a(1) = c(1) * x(1);
    a(2) = c(2) * x(1) * (x(1)-1) / 2;
    a(3) = c(3) * x(2);
    a(4) = c(4) * x(2);
end

function J = compute_jacobian(x, c, V, tau)
    % Jacobian of G(X) = X - V*(a(X)*tau) - RHS
    J = eye(3); % Identity matrix for dX/dX
    % Derivatives of propensities
    da1 = [c(1), 0, 0];
    da2 = [c(2)*(2*x(1)-1)/2, 0, 0];
    da3 = [0, c(3), 0];
    da4 = [0, c(4), 0];
    % Contribution from each reaction
    J = J - tau * [...
        V(:,1)*da1 + V(:,2)*da2 + V(:,3)*da3 + V(:,4)*da4...
    ];
end

function [tau, critical] = select_tau(x, c, epsilon, delta, remaining_time)
    % Simplified tau selection focusing on non-equilibrium reactions
    %a = compute_propensities(x, c);
    g = [2 + 1/(x(1)-1), 1, 1]; % Example order parameters
    
    % For R1 (S1 consumption)
    mu1 = -c(1)*x(1);
    var1 = c(1)*x(1);
    tau1 = min(...
        max(epsilon*x(1)/g(1), 1)/abs(mu1), ...
        (max(epsilon*x(1)/g(1), 1))^2 / var1 ...
    );
    
    % For R4 (S2 consumption)
    mu4 = -c(4)*x(2);
    var4 = c(4)*x(2);
    tau4 = min(...
        max(epsilon*x(2)/g(2), 1)/abs(mu4), ...
        (max(epsilon*x(2)/g(2), 1))^2 / var4 ...
    );
    
    tau = min([tau1, tau4, remaining_time]);
    critical = false;
end

function plot_results(results)
    figure('Color','white');
    plot(results(:,1), results(:,2), 'LineWidth', 1.5, 'Color', [0.12,0.47,0.71]);
    hold on;
    plot(results(:,1), results(:,3), 'LineWidth', 1.5, 'Color', [0.84,0.37,0.00]);
    plot(results(:,1), results(:,4), 'LineWidth', 1.5, 'Color', [0.17,0.63,0.17]);
    hold off;
    legend({'S_1', 'S_2', 'S_3'}, 'FontSize', 12, 'Location', 'best');
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Molecule Count', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 12);
end