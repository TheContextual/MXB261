% some cheeky parameters
exposedToInfectedRate = 1/5; % beta
transmissionRate = 10;       % alpha
recoveryRate = 1/10;         % rho

initialSusceptible = 990;    %S0
initialExposed = 0;          %E0
initialInfected = 10;        %I0
initialRecovered = 0;        %R0

% T = period
T = 360;

% h = step size for Euler's method
h = 0.5;

% N = population size
N = initialSusceptible + initialExposed + initialInfected + initialRecovered;

% this is the 'time vector' a vector of values between 0 and the period 
t = 0:h:T;
numSteps = length(t); % number of steps for Euler's method

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

% matlab is complaining about 't' but even though its not used, you kinda need it for ode calculations
function dydt = seir_ode(t, y, exposedToInfectRate, transmissionRate, recoveryRate) 

% y holds the information on each population in a vector as follows 
S = y(1); 
E = y(2); 
I = y(3); 
R = y(4); 

% this is total population as usual (assuming no births or deaths of course)
N = S + E + I + R; 

% same calcs as in Euler method, reasoning is the same
dSdt = -transmissionRate * S * I / N; 
dEdt = transmissionRate * S * I / N - exposedToInfectRate * E; 
dIdt = exposedToInfectRate * E - recoveryRate * I;
dRdt = recoveryRate * I;  

% at the end here we are just returning these derivatives as a column vector (thats just how ode works)  
dydt = [dSdt; dEdt; dIdt; dRdt]; 
end

% here we actually do the ode calculation
[t_ode, Y_ode] = ode45(@(t, y) seir_ode(t, y, exposedToInfectedRate, transmissionRate, recoveryRate), [0 T], [initialSusceptible, initialExposed, initialInfected, initialRecovered]);

% aaannddd we plot everything
figure;
hold on;
plot(t, S, 'g', t, E, 'b', t, I, 'r', t, R, 'k');
plot(t_ode, Y_ode(:, 1), 'g--', t_ode, Y_ode(:, 2), 'b--', t_ode, Y_ode(:, 3), 'r--', t_ode, Y_ode(:, 4), 'k--');
xlabel('Time (days)');
ylabel('Population');
legend('Euler Susceptible', 'Euler Exposed', 'Euler Infected', 'Euler Recovered', ...
       'ode45 Susceptible', 'ode45 Exposed', 'ode45 Infected', 'ode45 Recovered');
title('Comparison of Euler''s Method and ode45');
hold off;

% chosen accuracy measure will be RMSE

% we need to interpolate the Euler results according to the same points in time as ode45, because they need to be directly compared
S_interp = interp1(t, S, t_ode); 
E_interp = interp1(t, E, t_ode); 
I_interp = interp1(t, I, t_ode); 
R_interp = interp1(t, R, t_ode); 

% RMSE for each
RMSE_S = sqrt(mean((S_interp - Y_ode(:, 1)).^2))
RMSE_E = sqrt(mean((E_interp - Y_ode(:, 2)).^2)) 
RMSE_I = sqrt(mean((I_interp - Y_ode(:, 3)).^2)) 
RMSE_R = sqrt(mean((R_interp - Y_ode(:, 4)).^2))