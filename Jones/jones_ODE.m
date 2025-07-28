function jones_ODE()
    % Parameters
    KK = 3.36;
    NC = 4.11;
    NN = 285.0;
    a = 1.6;
    alpha = 0.195;
    c = 13.0;
    d = 0.01;
    delta = 0.7;
    gamma = 1.7475e-6;
    k = 1.7437e-7;
    Lambda = 10000.0;
    mu = 0.07;

    % Initial conditions: [A, CC, T, Tstar, V]
    X0 = [1.0, 18800.0, 323000.0, 7800.0, 120000.0];
    tspan = linspace(0, 250, 1000);
    [t, X] = ode45(@(t, X) jones_odes(t, X, KK, NC, NN, a, alpha, c, d, delta, gamma, k, Lambda, mu), tspan, X0);
    species_names = {'A', 'CC', 'T', 'T^*', 'V'};
    figure('Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);
    hold on;
    colororder(lines(length(species_names)));  % MATLAB colormap

    for i = 1:length(species_names)
        plot(t, X(:, i), 'LineWidth', 2.5);
    end

    hold off;
    xlabel('Time (days)', 'FontSize', 16);
    ylabel('Value (Concentration)', 'FontSize', 16);
    legend(species_names, 'Location', 'northeastoutside', 'FontSize', 12);
    set(gca, 'FontSize', 14);
    grid on;
end

function dX = jones_odes(~, X, KK, NC, NN, a, alpha, c, d, delta, gamma, k, Lambda, mu)
    A = X(1);
    CC = X(2);
    T = X(3);
    Tstar = X(4);
    V = X(5);
    dA = -gamma*A*T;
    dCC = alpha*k*T*V - mu*CC;
    dT = Lambda + (a*A*T)/(KK + A) - k*T*V - d*T;
    dTstar = k*T*V - delta*Tstar - alpha*k*T*V;
    dV = delta*NN*Tstar + mu*NC*CC - c*V;

    dX = [dA; dCC; dT; dTstar; dV];
end
