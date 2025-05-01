function RRE_MM()
    % Constants
    vol = 1e-15;
    n_A = 6.023e23;
    k = [1e6 / (n_A * vol), 1e-4, 1e-1];  % [k1, k2, k3]
    initials = [round(5e-7 * n_A * vol), round(2e-7 * n_A * vol), 0, 0];
    t_end = 50;
    t = linspace(0, t_end, 1000);

    % Solve ODEs
    [t, y] = ode45(@(t,y) michaelis_menten_system(y, k), t, initials);
    results = [t, y];

    % Plot
    figure('Color', 'white');
    plot(results(:,1), results(:,2), 'LineWidth', 3, 'Color', '#1f77b4'); % S
    hold on;
    plot(results(:,1), results(:,4), 'LineWidth', 3, 'Color', '#2ca02c'); % ES
    plot(results(:,1), results(:,3), 'LineWidth', 3, 'Color', '#ff7f0e'); % E
    plot(results(:,1), results(:,5), 'LineWidth', 3, 'Color', '#d62728'); % P
    hold off;
    
    legend({'Substrate (S)', 'Complex (ES)', 'Enzyme (E)', 'Product (P)'}, ...
           'Location', 'best', 'FontSize', 16);
    xlabel('Time (s)', 'FontSize', 14); 
    ylabel('Molecule count', 'FontSize', 14);
    grid on;
    set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
end

function dydt = michaelis_menten_system(y, k)
    % Unpack variables
    s = y(1);  % (S)
    e = y(2);  % (E)
    c = y(3);  % (ES)
    %derivatives
    ds_dt = -k(1) * e * s + k(2) * c;
    de_dt = -k(1) * e * s + (k(2) + k(3)) * c;
    dc_dt =  k(1) * e * s - (k(2) + k(3)) * c;
    dp_dt =  k(3) * c;
    dydt = [ds_dt; de_dt; dc_dt; dp_dt];
end
