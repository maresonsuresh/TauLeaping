function DDM_SSA_stars()
    % Rate constants
    C1 = 0.05;     % R1: S1 → ∅
    C2 = 50;       % R2: 2S1 → S2
    C3 = 1e6;      % R3: S2 → 2S1
    C4 = 0.05;     % R4: S2 → S3
    c = [C1, C2, C3, C4];

    % Stoichiometric matrix (4 reactions, 3 species)
    V = [-1  -2  +2   0;   % ΔS1
          0  +1  -1  -1;   % ΔS2
          0   0   0  +1];  % ΔS3

    % Initial state [S1, S2, S3]
    X = [400; 800; 0];
    t = 0;
    tfinal = 1;

    % Preallocate result storage
    max_steps = 3*1e8;
    results = NaN(max_steps, 4); % Columns: [time, S1, S2, S3]
    ii = 1;
    results(ii,:) = [t X'];

    % Gillespie SSA
    rng('shuffle');
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

    % Uniform time points for interpolation
    n_points = 30;  % you can change this for denser/sparser stars
    t_uniform = linspace(results(1,1), results(end,1), n_points);
    s1_interp = interp1(results(:,1), results(:,2), t_uniform, 'previous');
    s2_interp = interp1(results(:,1), results(:,3), t_uniform, 'previous');

    % Plot interpolated values with filled bold stars and connecting lines
    figure('Color','white');
    plot(t_uniform, s1_interp, '-p',  'LineWidth', 1.5,...
       'MarkerSize', 8, 'Color', '#1f77b4', ...
       'MarkerFaceColor', '#1f77b4');  % S1

    hold on;
    plot(t_uniform, s2_interp, '-p', 'LineWidth', 1.5, ...
        'MarkerSize', 8, 'Color', '#ff7f0e', ...
        'MarkerFaceColor', '#ff7f0e');  % S2

    % Set the x-axis limits to make sure the plot stops exactly at t = 1
    xlim([0, tfinal]);

    legend({'S_1', 'S_2'}, 'FontSize', 12, 'Location', 'best');
    xlabel('Time (s)');
    ylabel('Molecule count');
    grid on;
    set(gca, 'FontSize', 12);
end
