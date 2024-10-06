%% Initial parameters
exposedToInfectedRate = 1/5; % beta
transmissionRate = 1/14;     % alpha
recoveryRate = 1/10;         % rho

S0 = 990;    % Initial susceptible individual
E0 = 0;      % Initial exposed individual
I0 = 10;     % Initial infected individual
R0 = 0;      % Initial recovered individual
Population = S0 + E0 + I0 + R0;

% T = period
T = 360;

% State vector
x = [S0;
    E0;
    I0;
    R0];

% initialise time and state vector for tracking
t(1) = 0;
state_vec_tracking(1,:) = x;

% Run simulations
for i = 1:t
    % unfinished bwaa
end
