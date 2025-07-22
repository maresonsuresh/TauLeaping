function Collins_IPTG_ATC()
% Parameters
alpha1 = 28.98;
alpha2 = 28.98;
k = 300;
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

% Preallocate results
max_steps = 3.5e8;
Results = zeros(max_steps, 5);
Results(1, :) = [t, X'];

step = 1;

tic;
while t < Tmax
    p1 = X(1); p2 = X(2);
    m1 = X(3); m2 = X(4);

    beta1 = 1;
    beta2 = 1;
    if t < 500
        beta1 = 1;
        beta2 = 0.3;
    elseif t < 1000
        beta1 = 0.3;
        beta2 = 1;
    elseif t < 3000
        beta1 = 1;
        beta2 = 0.3;
    elseif t < 5000
        beta1 = 0.3;
        beta2 = 1;
    elseif t < 7500
        beta1 = 1;
        beta2 = 0.3;
    elseif t < 7900
        beta1 = 0.3;
        beta2 = 1;
    end
       

    % Propensities
    a(1) = k * beta1 * alpha1;
    a(2) = k * C * m1;
    a(3) = C * m1;
    a(4) = C * p1;
    a(5) = k * beta2 * alpha2;
    a(6) = k * C * m2;
    a(7) = C * m2;
    a(8) = C * p2;

    a0 = sum(a);
    if a0 == 0, break; end

    tau = -log(rand) / a0;
    t = t + tau;

    r = rand * a0;
    j = find(cumsum(a) >= r, 1);
    X = X + V(:,j);
    X = max(X, 0); % non-negative

    step = step + 1;
    Results(step, :) = [t, X'];
end

% Trim results
Results = Results(1:step, :);
time = Results(:,1);
p1 = Results(:,2);
p2 = Results(:,3);
m1 = Results(:,4);
m2 = Results(:,5);

elapsed_time = toc;
fprintf('Simulation time: %.4f seconds\n', elapsed_time);
fprintf('Number of recorded steps: %d\n', size(Results, 1));

% Condensed plot using tiledlayout
figure;
tiledlayout(2,2,'Padding','tight','TileSpacing','tight');
data = [p1, p2, m1, m2];
titles = {'Protein 1','Protein 2','mRNA 1','mRNA 2'};
for i = 1:4
    nexttile; plot(time, data(:,i));
    title(titles{i}); xlim([0 Tmax]); ylim([0 200]);
end
end
