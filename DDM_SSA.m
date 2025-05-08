function DDM_SSA()
    % Rate constants + Stoichiometric matrix
    c = [1, 10, 1000, 0.1];
    V = [-1  -2  +2   0;   % S1
          0  +1  -1  -1;   % S2
          0   0   0  +1];  % S3
    % Initial state [S1, S2, S3]
    X = [400; 798; 0];
    t = 0;
    tfinal = 0.2;

    % Preallocate result storage
    max_steps = 5*1e8;
    results = NaN(max_steps, 4); % Columns: [time, S1, S2, S3]
    ii = 1;
    results(ii,:) = [t X'];

    % Gillespie SSA
    rng(45);
    tic;
    while t < tfinal
        % Propensities
        a = [c(1)*X(1);                     % R1: S1 → ∅
             c(2)*X(1)*(X(1)-1)/2;          % R2: 2S1 → S2
             c(3)*X(2);                     % R3: S2 → 2S1
             c(4)*X(2)];                    % R4: S2 → S3
        a0 = sum(a);

        % Stop if no more reactions can occur
        if a0 == 0
            warning('All propensities zero. Simulation halted at time %.5f.', t);
            break;
        end

        % Time to next reaction
        tau = log(1/rand)/a0;
        % Select reaction
        j = find(cumsum(a) >= rand*a0, 1);

        % Update time and state
        t = t + tau;
        X = X + V(:,j);

        ii = ii + 1;
        if ii > max_steps
            warning('Maximum number of steps (%d) reached. Ending simulation early.', max_steps);
            break;
        end

        results(ii,:) = [t X'];
    end
    elapsed_time = toc;
    fprintf('Simulation completed in %.4f seconds with %d steps.\n', elapsed_time, ii);

    % Trim unused rows
    results = results(1:ii, :);
    % Plotting
    figure('Color','white');
    plot(results(:,1), results(:,2), 'LineWidth', 1, 'Color', '#1f77b4'); % S1
    hold on;
    plot(results(:,1), results(:,3), 'LineWidth', 1, 'Color', '#ff7f0e'); % S2
    plot(results(:,1), results(:,4), 'LineWidth', 1, 'Color', '#2ca02c'); % S3
    hold off;
    legend({'S_1', 'S_2', 'S_3'}, 'FontSize', 12, 'Location', 'best');
    xlabel('Time (s)');
    ylabel('Molecule count');
    grid on;
    set(gca, 'FontSize', 12);
end
