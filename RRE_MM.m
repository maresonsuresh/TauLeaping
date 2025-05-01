function RRE_MM()
    % Constants
    vol = 1e-15;
    n_A = 6.023e23;
    c1 = 1e6 / (n_A * vol);
    c2 = 1e-4;
    c3 = 1e-1;
    S_initial = round(5e-7 * n_A * vol);
    E_initial = round(2e-7 * n_A * vol);
    ES_initial = 0;
    P_initial = 0;
    t_end = 50;
    t = linspace(0, t_end, 1000);

    % Solve ODEs
    [t, y] = ode45(@(t,y) michaelis_menten_system(t, y, c1, c2, c3, E_initial), t, [S_initial; ES_initial; E_initial; P_initial]);
    results = [t, y(:,1), y(:,3), y(:,2), y(:,4)];

    % Plot
    figure('Color', 'white');
    plot(results(:,1), results(:,2), 'LineWidth', 3, 'Color', '#1f77b4'); % S
    hold on;
    plot(results(:,1), results(:,4), 'LineWidth', 3, 'Color', '#2ca02c'); % ES
    plot(results(:,1), results(:,3), 'LineWidth', 3, 'Color', '#ff7f0e'); % E
    plot(results(:,1), results(:,5), 'LineWidth', 3, 'Color', '#d62728'); % P
    hold off;
    
    legend({'Substrate (S)', 'Complex (ES)', 'Enzyme (E)', 'Product (P)'}, 'Location', 'best', 'FontSize', 16);
    xlabel('Time (s)', 'FontSize', 14); 
    ylabel('Molecule count', 'FontSize', 14);
    grid on;
    set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
end

function dydt = michaelis_menten_system(~, y, c1, c2, c3, E)
    s = y(1);
    c = y(2);
    e = y(3);
    p = y(4);
    
    ds_dt = -c1 * E * s + (c1 * s + c2) * c;
    dc_dt = c1 * E * s - (c1 * s + c2 + c3) * c;
    de_dt = -c1 * E * s + (c1 * s + c2) * c + c3 * c;
    dp_dt = c3 * c;
    dydt = [ds_dt; dc_dt; de_dt; dp_dt];
end
