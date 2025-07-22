function collins_memory()
    % Parameters
    alpha1 = 28.98;
    alpha2 = 28.98;
    k = 1000;
    C = 0.23;
    V = [0  0  1 -1  0  0  0  0;   % protein
         0  0  0  0  0  0  1 -1;   % protein
         1 -1  0  0  0  0  0  0;   % mRNA
         0  0  0  0  1 -1  0  0];  % mRNA
    rng(3);

    % Initial state: [p1, p2, m1, m2]
    X = [76; 75; 60; 60];

    % Time
    t = 0;
    Tmax = 8000;

    % Memory-efficient: sample only at integer times
    num_points = Tmax + 1;
    sampled_times = (0:Tmax)';
    sampled_results = zeros(num_points, 5);  % [t, p1, p2, m1, m2]
    sampled_results(1, :) = [t, X'];
    next_sample_index = 2;

    tic;
    while t < Tmax
        p1 = X(1); p2 = X(2);
        m1 = X(3); m2 = X(4);

        % Time-dependent input
        beta1 = 1; beta2 = 1;
        if t < 500
            beta2 = 0.3;
        elseif t < 1000
            beta1 = 0.3; beta2 = 1;
        elseif t < 3000
            beta1 = 1; beta2 = 0.3;
        elseif t < 5000
            beta1 = 0.3; beta2 = 1;
        elseif t < 7500
            beta1 = 1; beta2 = 0.3;
        elseif t < 7900
            beta1 = 0.3; beta2 = 1;
        end

        % Propensities
        a = [
            k * beta1 * alpha1;
            k * C * m1;
            C * m1;
            C * p1;
            k * beta2 * alpha2;
            k * C * m2;
            C * m2;
            C * p2
        ];
        a0 = sum(a);
        if a0 == 0, break; end

        % Gillespie step
        tau = -log(rand) / a0;
        t_next = t + tau;
        r = rand * a0;
        j = find(cumsum(a) >= r, 1);

        % Record all integer times in [t, t_next)
        while next_sample_index <= num_points && sampled_times(next_sample_index) < t_next
            sampled_results(next_sample_index, :) = [sampled_times(next_sample_index), X'];
            next_sample_index = next_sample_index + 1;
        end

        % Update state
        X = X + V(:, j);
        X = max(X, 0);
        t = t_next;
    end

    % Fill in any remaining times with last state
    while next_sample_index <= num_points
        sampled_results(next_sample_index, :) = [sampled_times(next_sample_index), X'];
        next_sample_index = next_sample_index + 1;
    end

    elapsed_time = toc;
    fprintf('Simulation time: %.4f seconds\n', elapsed_time);
    fprintf('Final time reached: %.2f\n', t);

    % Plot
    time = sampled_results(:, 1);
    p1 = sampled_results(:, 2);
    p2 = sampled_results(:, 3);
    m1 = sampled_results(:, 4);
    m2 = sampled_results(:, 5);

    figure;
    tiledlayout(2,2,'Padding','tight','TileSpacing','tight');
    data = [p1, p2, m1, m2];
    titles = {'Protein 1','Protein 2','mRNA 1','mRNA 2'};
    for i = 1:4
        nexttile; plot(time, data(:,i));
        title(titles{i}); xlim([0 Tmax]); ylim([0 200]);
    end
end
