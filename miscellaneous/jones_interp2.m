function jones()
    format long g
    % Parameters
    gamma = 1.7475e-06;
    lambda = 10000;
    a = 1.6;
    KK = 3.36;
    k = 1.7437e-07;
    d = 0.01;
    delta = 0.7;
    mu = 0.07;
    alpha = 0.195;
    c = 13;
    NN = 285;
    NC = 4.11;

    % Initial state: [A, CC, T, Tstar, V]
    X = [1; 18800; 323000; 7800; 120000]; 
    % Stoich
    Q = [
    -1   0   0   0   0   0   0   0   0   0   0;  % A
     0   0   0   0   0   0   0   0  -1   0   1;  % CC
     0   0   0   1   1  -1  -1   0   0   0   0;  % T
     0   0   0   0   0   1   0  -1   0   0  -1;  % Tstar
     0   1   1   0   0   0   0   0   0  -1   0   % V
    ];

    % Settings
    rng(1);
    t = 0;
    tfinal = 250;
    max_steps = 1e9;
    
    % Preallocate for 251 points: t=0,1,2,...,250
    num_points = 251;
    sampled_times = (0:250)';  % Column vector [0;1;...;250]
    sampled_results = zeros(num_points, 6); % [t, A, CC, T, Tstar, V]
    next_sample_index = 1;  % Start with first sample (t=0)
    
    % Record initial state (t=0)
    sampled_results(next_sample_index, :) = [t, X'];
    next_sample_index = 2;  % Next sample at t=1
    
    step = 0;
    toc_start = tic;
    
    % Gillespie simulation
    while t < tfinal && step < max_steps
        A = X(1); CC = X(2); T = X(3); Tstar = X(4); V = X(5);
        
        % Reaction propensities
        a_vec = [
            gamma*A*T;           % v1
            delta*NN*Tstar;      % v10
            mu*NC*CC;            % v11
            lambda;              % v2
            a*A*T/(KK + A);      % v3
            k*T*V;               % v4
            d*T;                 % v5
            delta*Tstar;         % v6
            mu*CC;               % v7
            c*V;                 % v8
            alpha*k*T*V;         % v9
        ];
        
        a0 = sum(a_vec);
        if a0 == 0
            break;  % Exit if no reactions can occur
        end
        
        % Time to next event and reaction selection
        tau = -log(rand) / a0;
        r = rand * a0;
        cum_a = cumsum(a_vec);
        j = find(r < cum_a, 1, 'first');
        t_next = t + tau;
        
        % Record states at integer times in [t, t_next)
        while next_sample_index <= num_points
            t_sample = sampled_times(next_sample_index);
            if t_sample < t_next
                % State constant in [t, t_next) → use current state
                sampled_results(next_sample_index, :) = [t_sample, X'];
                next_sample_index = next_sample_index + 1;
            else
                break;
            end
        end
        
        % Stop if next event is beyond tfinal
        if t_next >= tfinal
            break;
        end
        
        % Update state for the reaction
        X = X + Q(:, j);
        X = max(X, 0);  % Ensure non-negativity
        t = t_next;
        step = step + 1;
    end
    
    % Fill remaining sample points with last state
    while next_sample_index <= num_points
        t_sample = sampled_times(next_sample_index);
        sampled_results(next_sample_index, :) = [t_sample, X'];
        next_sample_index = next_sample_index + 1;
    end
    
    % Output performance metrics
    elapsed_time = toc(toc_start);
    fprintf('Simulation time: %.4f seconds\n', elapsed_time);
    fprintf('Number of steps (events): %d\n', step);
    disp('Final state at t=250:');
    disp(sampled_results(end, :));
    
    % Plot results using the 251 points
    figure;
    hold on;
    plot(sampled_results(:,1), sampled_results(:,2), 'b-', 'LineWidth', 2, 'DisplayName', 'A');
    plot(sampled_results(:,1), sampled_results(:,3), 'y-', 'LineWidth', 2, 'DisplayName', 'CC');
    plot(sampled_results(:,1), sampled_results(:,4), 'g-', 'LineWidth', 2, 'DisplayName', 'T');
    plot(sampled_results(:,1), sampled_results(:,5), 'r-', 'LineWidth', 2, 'DisplayName', 'T^*');
    plot(sampled_results(:,1), sampled_results(:,6), 'm-', 'LineWidth', 2, 'DisplayName', 'V');
    hold off;
    legend('Location', 'best');
    ylabel('Molecule Count');
    grid on;
end