function SSA_MM()
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
    tfinal = 50;
    results = [t X'];
    
    % Gillespie
    while t < tfinal 
        % Step 1: Evaluate propensities and their sum
        a = [c(1)*X(1)*X(2)   % S+E->ES
             c(2)*X(3)         % ES->S+E
             c(3)*X(3)];      % ES->P+E
        a0 = sum(a);

        % Step 3: Select reaction index j
        j = find(cumsum(a) >= rand*a0, 1);
        % Step 4: Time step
        tau = log(1/rand)/a0;
        % Step 5: Update state, time, and results
        X = X + V(:,j);
        t = t + tau;
        results(end+1,:) = [t X'];
    end
    
    % Plot
    figure('Color','white');
    plot(results(:,1), results(:,2), 'LineWidth', 2, 'Color', '#1f77b4'); % S
    hold on;
    plot(results(:,1), results(:,4), 'LineWidth', 2, 'Color', '#2ca02c'); % ES
    plot(results(:,1), results(:,3), 'LineWidth', 2, 'Color', '#ff7f0e'); % E
    plot(results(:,1), results(:,5), 'LineWidth', 2, 'Color', '#d62728'); % P
    hold off;
    
    legend({'Substrate (S)','Complex (ES)','Enzyme (E)','Product (P)'}, ...
           'Location', 'best');
    xlabel('Time (s)'); 
    ylabel('Molecule count');
    title('Michaelis-Menten SSA Simulation');
    grid on;
    set(gca, 'FontSize', 12);
end