function ToggleSwitch_SSA()
    % Parameters from the table
    params = struct(...
        'km0L', 0.3045, ...
        'km0T', 0.3313, ...
        'kmL', 13.01, ...
        'kmT', 5.055, ...
        'kpL', 0.6606, ...
        'kpT', 0.5098, ...
        'gmL', 0.1386, ...
        'gmT', 0.1386, ...
        'gpL', 0.0165, ...
        'gpT', 0.0165, ...
        'thetaLacI', 124.9, ...
        'thetaTetR', 76.40, ...
        'etaLacI', 2.00, ...
        'etaTetR', 2.152 ...
    );

    % Initial state [LacI_protein, TetR_protein, LacI_mRNA, TetR_mRNA]
    X = [10; 10; 0; 0];
    t = 0;
    tfinal = 3999;
    max_steps = 1e6;
    results = zeros(max_steps, 5); % [time, LacI_p, TetR_p, mRNA_L, mRNA_T]
    results(1,:) = [t, X'];
    step = 1;

    % Stoichiometry matrix
    V = [0  0  1 -1  0  0  0  0;   % LacI protein
         0  0  0  0  0  0  1 -1;   % TetR protein
         1 -1  0  0  0  0  0  0;   % LacI mRNA
         0  0  0  0  1 -1  0  0];  % TetR mRNA

    rng('shuffle');

    while t < tfinal && step < max_steps
        LacIp = X(1); TetRp = X(2);
        mLacI = X(3); mTetR = X(4);

        % Propensities
        a = [params.km0L + params.kmL / (1 + (TetRp / params.thetaTetR)^params.etaTetR); % LacI transcription
             params.gmL * mLacI;                                                         % LacI mRNA decay
             params.kpL * mLacI;                                                         % LacI translation
             params.gpL * LacIp;                                                         % LacI protein decay
             params.km0T + params.kmT / (1 + (LacIp / params.thetaLacI)^params.etaLacI); % TetR transcription
             params.gmT * mTetR;                                                         % TetR mRNA decay
             params.kpT * mTetR;                                                         % TetR translation
             params.gpT * TetRp];                                                        % TetR protein decay

        a0 = sum(a);
        if a0 == 0
            break;
        end

        tau = -log(rand) / a0;
        r = rand * a0;
        j = find(cumsum(a) >= r, 1);

        % Update time and state
        t = t + tau;
        X = X + V(:,j);
        step = step + 1;
        results(step, :) = [t, X'];
    end

    results = results(1:step, :);

    % Plot
    figure;
yyaxis left
plot(results(:,1), results(:,2), 'b', 'LineWidth', 1); % LacI
ylabel('LacI Protein (blue)');
ylim([0, 3300]);

yyaxis right
plot(results(:,1), results(:,3), 'r', 'LineWidth', 1); % TetR
ylabel('TetR Protein (red)');
ylim([0, 220]);

xlabel('Time (min)');
title('LacI and TetR Protein Dynamics');
grid on;
end
