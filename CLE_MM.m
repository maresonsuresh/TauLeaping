function CLE_MM()
    % Parameters and Initial Conditions
    nA = 6.023e23;
    vol = 1e-15;
    c = [1e6/(nA*vol); 1e-4; 0.1];
    % Stoichiometric matrix
    V = [-1  1  0; -1  1  1; 1 -1 -1; 0  0  1];
    % Initial state
    X = [round(5e-7*nA*vol); round(2e-7*nA*vol); 0; 0];
    % Simulation setup
    t = 0;
    tfinal = 49;
    results = [t X'];
    % Fixed time step for tau-leaping 
    tau = 1; 
    
    tic;
    % CLE simulation
    while t < tfinal
        
        % Evaluate propensities
        a = [c(1)*X(1)*X(2)   % S+E->ES
             c(2)*X(3)         % ES->S+E
             c(3)*X(3)];       % ES->P+E
        
        % Generate normal random variables / update state
        xi = randn(3,1);
        X_new = X + V(:,1:3) * a * tau + V(:,1:3) * (sqrt(a .* tau) .* xi);
        
        % Update state, time, and results
        X = X_new;
        t = t + tau;
        results(end+1,:) = [t X'];
    end
    
    elapsed_time = toc;
    fprintf('Simulation time: %.4f seconds\n', elapsed_time);

    % Plot
    figure('Color','white');
    plot(results(:,1), results(:,2), 'LineWidth', 3, 'Color', '#1f77b4'); % S
    hold on;
    plot(results(:,1), results(:,4), 'LineWidth', 3, 'Color', '#2ca02c'); % ES
    plot(results(:,1), results(:,3), 'LineWidth', 3, 'Color', '#ff7f0e'); % E
    plot(results(:,1), results(:,5), 'LineWidth', 3, 'Color', '#d62728'); % P
    hold off;
    
    legend({'Substrate (S)','Complex (ES)','Enzyme (E)','Product (P)'}, ...
           'Location', 'best', 'FontSize', 16);
    xlabel('Time (s)'); 
    ylabel('Molecule count');
    grid on;
    set(gca, 'FontSize', 12);
end
