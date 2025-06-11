% Mechanisms of noise-resistance in genetic oscillators
% https://www.pnas.org/doi/epdf/10.1073/pnas.092133899
% https://www.ebi.ac.uk/biomodels/BIOMD0000000035

function gillespie_vilar_oscillator()
    % Reaction rates
    alphaA = 50;
    alphaAp = 500;
    alphaR = 0.01;
    alphaRp = 50;
    betaA = 50;
    betaR = 5;
    gammaA = 1;
    gammaC = 2;
    gammaR = 1;
    deltaA = 1;
    deltaMA = 10;
    deltaMR = 0.5;
    deltaR = 0.2;
    thetaA = 50;
    thetaR = 100;

    % Initial conditions
    initial_state = [0; 0; 1; 0; 1; 0; 0; 0; 0]; % [A; C; DA; DAp; DR; DRp; MA; MR; R]
    
    % Species names for plotting
    species_names = {'A', 'C', 'DA', 'DAp', 'DR', 'DRp', 'MA', 'MR', 'R'};
    
    % Stoichiometry matrix (columns = reactions, rows = species)
    nu = [0  0  0  0  1  0 -1 -1 -1 -1  0  0  0  0  1  1;  % A
          0  0  0  0  0  0  0  1  0  0 -1  0  0  0  0  0;  % C
          0  0  0  0  0  0 -1  0  0  0  0  0  0  0  1  0;  % DA
          0  0  0  0  0  0  1  0  0  0  0  0  0  0 -1  0;  % DAp
          0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  1;  % DR
          0  0  0  0  0  0  0  0  1  0  0  0  0  0  0 -1;  % DRp
          1  1  0  0  0  0  0  0  0  0  0 -1  0  0  0  0;  % MA
          0  0  1  1  0  0  0  0  0  0  0  0 -1  0  0  0;  % MR
          0  0  0  0  0  1  0 -1  0  0  1  0  0 -1  0  0]; % R
    
    % Simulation parameters
    t_max = 200; % Maximum simulation time
    
    % Run Gillespie simulation (records every state change)
    [time, states] = gillespie_simulation_continuous(initial_state, nu, t_max, ...
        alphaA, alphaAp, alphaR, alphaRp, betaA, betaR, gammaA, gammaC, gammaR, ...
        deltaA, deltaMA, deltaMR, deltaR, thetaA, thetaR);
    
    % Plot results
    plot_all_in_one(time, states, species_names);
end

function [time, states] = gillespie_simulation_continuous(initial_state, nu, t_max, ...
    alphaA, alphaAp, alphaR, alphaRp, betaA, betaR, gammaA, gammaC, gammaR, ...
    deltaA, deltaMA, deltaMR, deltaR, thetaA, thetaR)
    
    % Initialize
    state = initial_state;
    time = 0;
    
    % Preallocate with reasonable size (will grow if needed)
    states = zeros(1000000, length(initial_state));
    times = zeros(1000000, 1);
    states(1,:) = state';
    times(1) = time;
    record_index = 2;
    
    % Main simulation loop
    while time < t_max
        % Calculate reaction propensities
        A = state(1); C = state(2); DA = state(3); DAp = state(4);
        DR = state(5); DRp = state(6); MA = state(7); MR = state(8); R = state(9);
        
        a = [alphaA * DA;
             alphaAp * DAp;
             alphaR * DR;
             alphaRp * DRp;
             betaA * MA;
             betaR * MR;
             gammaA * A * DA;
             gammaC * A * R;
             gammaR * A * DR;
             deltaA * A;
             deltaA * C;
             deltaMA * MA;
             deltaMR * MR;
             deltaR * R;
             thetaA * DAp;
             thetaR * DRp];
        
        a0 = sum(a);
        
        if a0 == 0
            break; % No more reactions can occur
        end
        
        % Determine time to next reaction
        tau = -log(rand) / a0;
        time = time + tau;
        
        % If we've exceeded t_max, adjust tau to end exactly at t_max
        if time > t_max
            tau = tau - (time - t_max);
            time = t_max;
        end
        
        % Determine which reaction occurs
        r = find(cumsum(a) >= rand * a0, 1);
        
        % Update state
        state = state + nu(:, r);
        
        % Ensure no negative values (can happen due to discrete nature)
        state(state < 0) = 0;
        
        % Record the new state
        states(record_index, :) = state';
        times(record_index) = time;
        record_index = record_index + 1;
        
        if time >= t_max
            break;
        end
    end
    
    % Trim unused preallocated space
    states = states(1:record_index-1, :);
    time = times(1:record_index-1);
end

function plot_all_in_one(time, states, species_names)
    figure;
    hold on;
    
    % Define a color palette
    colors = lines(9); % Get 9 distinct colors
    
    % Plot all species with solid lines
    for i = 1:9
        plot(time, states(:, i), 'LineWidth', 1.5, 'DisplayName', species_names{i}, 'Color', colors(i,:));
    end
    
    hold off;
    xlabel('Time');
    ylabel('Molecules');
    legend('Location', 'best');
    grid on;
    box on;
end