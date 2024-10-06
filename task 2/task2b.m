%% Initial parameters
exposedToInfectedRate = 1/5; % beta
transmissionRate = 1/14;     % alpha
recoveryRate = 1/10;         % rho

S0 = 990;    % Initial susceptible individual
E0 = 0;      % Initial exposed individual
I0 = 10;     % Initial infected individual
R0 = 0;      % Initial recovered individual
Population = S0 + E0 + I0 + R0;

T = 360; % Period

t = 1; % initial time

x = [S0;
    E0;
    I0;
    R0]; % State vector

% Stoichiometric vectors
v1 = [-1;
    +1;
    0;
    0];

v2 = [0;
    -1;
    +1;
    0];

v3 = [0;
    0;
    -1;
    +1];

% initialise time and state vector for tracking
t_vector(1) = 0;
state_vec_tracking(1,:) = x;

% Run simulations
for i = 1:2 % number of simulations
    while t < T
        count = 1;
        % Calculate q rates
        q1 = (exposedToInfectedRate*x(1)*x(3))/Population;
        q2 = transmissionRate*x(2);
        q3 = recoveryRate*x(3);

        propensity_vec = [q1, q2, q3];
        a0 = sum(propensity_vec);

        cumulative_prob = cumsum(propensity_vec) / a0;

        % Calculate wait time to next event
        u1 = rand;
        dt = -log(u1) / a0;
        % update time
        t = t + dt;
        t_vector(count + 1) = t;

        % determine which event occurred
        u2 = rand;

        if u2 < cumulative_prob(1) % Susceptible individual becomes exposed
            x = x + v1;
        elseif u2 < cumulative_prob(2) % Exposed individual is infected
            x = x + v2;
        elseif u2 < cumulative_prob(3) % Infected individual recovers
            x = x + v3;
        end
        state_vec_tracking(count+1,:) = x;
        count = count + 1;

        %plot
    end
end
