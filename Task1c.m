% some cheeky parameters
exposedToInfectedRate = 1/5; % beta
transmissionRate = 1/14;     % alpha
recoveryRate = 1/10;         % rho

initialSusceptible = 990;    %S0
initialExposed = 0;          %E0
initialInfected = 10;        %I0
initialRecovered = 0;        %R0

% T = period
T = 360;

% h = step size
h = 0.5;

% N = population size
N = initialSusceptible + initialExposed + initialInfected + initialRecovered;

% this is the 'time vector' a vector of values between 0 and the period 
t = 0:h:T;
numSteps = length(t);

% need to initialise our variables ay?
% S = Susceptible
% E = Exposed
% I = Infected
% R = Recovered
S = zeros(1, numSteps);
E = zeros(1, numSteps);
I = zeros(1, numSteps);
R = zeros(1, numSteps);

% setting those initial values to our initial variables 
S(1) = initialSusceptible;
E(1) = initialExposed;
I(1) = initialInfected;
R(1) = initialRecovered;

% Euler method
% we repeating this a select number of times
for k = 1:numSteps-1
    % the susceptible population should decrease with more contact with infected individuals
    S(k+1) = S(k) - h * transmissionRate * S(k) * I(k) / N;

    % the exposed population increases when susceptible become infected however decrease when they then become infectious themselves
    E(k+1) = E(k) + h * (transmissionRate * S(k) * I(k) / N - exposedToInfectedRate * E(k));

    % infected population simply increases when exposed become infectious and decrease when individuals recover
    I(k+1) = I(k) + h * (exposedToInfectedRate * E(k) - recoveryRate * I(k));

    % these recovered peoples then very simply are a result of infectious individuals completing the disease cycle
    R(k+1) = R(k) + h * (recoveryRate * I(k));
end

% aaannddd we plot everything
figure;
plot(t, S, 'g', t, E, 'b', t, I, 'r', t, R, 'k');
xlabel('Time (days)');
ylabel('Population');
legend('Susceptible', 'Exposed', 'Infected', 'Recovered');
title('SEIR Model Simulation using Euler''s Method');

% here is the equilibrium from part b so we can compare
S_eq = N 