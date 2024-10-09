% Code written by Rhys Lodding
% STUDENT NUMBER: n10757970
% MXB261 Group 11

clc; % Clear command window.

% Construct a grid of cells 101x101.
% Grid is split from (51, 1) to (51, 101) set as river tiles.
% Any agent cannot be on any river tile.
% River has [B] number of bridges that are placed randomly (treated as land
% tiles).
% Populate grid with 990 susceptible agents and 10 infected agents.
% Susceptible agents are randomly placed and infected agents are north side
% only.
% Only 1 agent per cell.
% Rules for Simulation:
% a) Each agent in turn will try to move up, down, left, or right.
%   - If agent goes out of bounds, no move happens.
%   - If agent goes on river tile, no move happens.
%   - If agent goes onto tile with another agent, no move happens.
%   - Agent can only move onto an empty land tile.
% b) After each agent moves, each infected agent randomly selects a
%    neighbouring cell. If cell contains susceptible agent then it becomes
%    exposed at a [alpha] percent chance.
% c) Each exposed agent will become an infected agent at a [beta] chance on
%    this step.
% d) Each infected agent will become recovered agent at a [rho] chance on
%    this step. (Recovered agents do not become infected again)
%
% Do simulations for two regimes and every multiple of 10 bridges up to
% 100.
% Regimes:
% [alpha] = 1, [beta] = 1/14, [rho] = 1/20
% [alpha] = 0.5, [beta] = 1/14, [rho] = 1/10
%
% Outcome:
% - Create a figure for each regime showing susceptible, exposed, infected, 
%   and recovered agents as a 2x5 subplot (for each number of bridges). 
%   Figure should go for 400 steps.
% - Create a movie for each plot no slower than 10 frames per second 
%   (40 seconds).

% Create Figure 1 for Regime 1
figure('Name', "Alpha = 1, Beta = 1/14, Rho = 1/20");
rows = 2;
cols = 5;
for index = 1:rows*cols
    subplot(rows, cols, index);
    runSimulationAndCreateSubplot(1, 1/14, 1/20, 10 * index, index, "Figure 1");
end

% Create Figure 2 for Regime 2
figure('Name', "Alpha = 0.5, Beta = 1/14, Rho = 1/10");
rows = 2;
cols = 5;
for index = 1:rows*cols
    subplot(rows, cols, index);
    runSimulationAndCreateSubplot(1, 1/14, 1/10, 10 * index, index, "Figure 2");
end

% Runs the spatial stochastic simulation
% alpha = chance susceptible agent gets turn exposed
% beta = chance exposed agents turn infected
% rho = chance infected agents turn recovered
% B = number of bridges across the river
% index = index of subplot for saving
% saveDir = directory of save file for images and movies
function runSimulationAndCreateSubplot(alpha, beta, rho, B, index, saveDir)
    fprintf("=== Running Simulation %d ===\n", index); % Debug

    % Create folders for images and videos
    saveFileDir = saveDir + '\\' + index;
    mkdir(saveFileDir);

    % Few Resused Settings
    gridSize = 101;
    numOfSusceptibleAgents = 990;
    numOfInfectedAgents = 10;
    numOfSteps = 400;
    frameMultiplier = 4; % Used to upscale the images and videos

    % Data of each type of agent for plots
    susceptibleAgentData = zeros(numOfSteps, 1);
    exposedAgentData = zeros(numOfSteps, 1);
    infectedAgentData = zeros(numOfSteps, 1);
    recoveredAgentData = zeros(numOfSteps, 1);

    % Start the video of the simulation
    video = createVideo(10, index, saveFileDir);
    
    % Create and populate city grid tiles
    cityGrid = createGrid(gridSize);
    cityGrid = addBridgesToGrid(cityGrid, B);
    fprintf("Created City Grid\n"); % Debug
    
    % Create and populate agent grid tiles
    agentGrid = createGrid(gridSize);
    agentGrid = populateGridWithAgents(agentGrid, cityGrid, numOfSusceptibleAgents, numOfInfectedAgents);
    fprintf("Populated Agents\n"); % Debug

    % Write first frame of video
    open(video);
    colourmap = createFrame(cityGrid, agentGrid, frameMultiplier);
    img_title = sprintf(saveFileDir + "\\start-%d.jpeg", index);
    imwrite(colourmap, img_title); % Debug
    writeVideo(video, colourmap); 

    printAgentDistribution(agentGrid); % Debug

    % Record number of each type of agent
    susceptibleAgentData(1) = countAllAgents(agentGrid, 1);
    exposedAgentData(1) = countAllAgents(agentGrid, 2);
    infectedAgentData(1) = countAllAgents(agentGrid, 3);
    recoveredAgentData(1) = countAllAgents(agentGrid, 4);
    
    % Loop through each step of simulation
    step = 1;
    while step <= numOfSteps
        agentGrid = randomlyMoveAgents(agentGrid, cityGrid); % Randomly move agents
        agentGrid = updateAgentStates(agentGrid, alpha, beta, rho); % Randomly update agents

        % Create frame for video
        colourmap = createFrame(cityGrid, agentGrid, frameMultiplier);
        writeVideo(video, colourmap);

        fprintf("Step %d Completed...\n", step); % Debug

        step = step + 1; % Iterate step
        printAgentDistribution(agentGrid); % Debug

        % Record number of each type of agent per step
        susceptibleAgentData(step) = countAllAgents(agentGrid, 1);
        exposedAgentData(step) = countAllAgents(agentGrid, 2);
        infectedAgentData(step) = countAllAgents(agentGrid, 3);
        recoveredAgentData(step) = countAllAgents(agentGrid, 4);
        
        % Stop simulation if no more exposed or infected
        if countAllAgents(agentGrid, 2) + countAllAgents(agentGrid, 3) == 0
            fprintf("Simulation Stopped As No More Infected Agents...\n"); % Debug

            % Cut off data to shorten plots
            susceptibleAgentData = susceptibleAgentData(1:step);
            exposedAgentData = exposedAgentData(1:step);
            infectedAgentData = infectedAgentData(1:step);
            recoveredAgentData = recoveredAgentData(1:step);

            break;
        end
    end

    colourmap = createFrame(cityGrid, agentGrid, frameMultiplier); % Debug
    img_title = sprintf(saveFileDir + "\\end-%d.jpeg", index);
    imwrite(colourmap, img_title); % Debug
    
    % Create the plot with data
    createPlot(susceptibleAgentData, exposedAgentData, infectedAgentData, recoveredAgentData, B);
    
    % Finish video
    close(video)
    fprintf("Simulation Complete!\n"); % Debug
end

% Creates a useable grid
% size = size of the grid
function grid = createGrid(size)
    grid = cell(size); % Create empty grid
    % Loop through all cells
    for x = 1:size
        for y = 1:size
            grid{x, y} = 0; % Set all cells to state 0
        end
    end
end

% Adds a river and bridges to grid
% grid = referenced grid
% numOfBridges = num of bridges across the river
function gridAfterBridges = addBridgesToGrid(grid, numOfBridges)
    % Enum for agents:
    % -1 = river cell
    % 0 = city cell
    % 1 = bridge cell

    bridgesCounter = 0;
    y_coord = round(length(grid) / 2); % Middle of the grid
    % Loop through each cell on x axis to set as river
    for x_coord = 1:length(grid)
        grid{y_coord, x_coord} = -1; % Set cell to water state (-1)
    end
    % Loop through each bridge to place it
    while (bridgesCounter < numOfBridges)
        x_coord = randi(length(grid)); % Get random x coordinate
        if grid{y_coord, x_coord} == -1 % Check if it's a valid river cell
            grid{y_coord, x_coord} = 1; % Set cell to bridge state (1)
            bridgesCounter = bridgesCounter + 1;
        end
    end
    gridAfterBridges = grid;
end

% Populated grid with agents
% grid = referenced grid
% cityGrid = grid of city tiles to condition spawning
% numOfSusceptible = number of agents susceptible to getting sick.
% numOfInfected = number of agents that can infect other agents.
function gridAfterPopulation = populateGridWithAgents(grid, cityGrid, numOfSusceptible, numOfInfected)
    % Enum for agents:
    % 0 = empty cell
    % 1 = susceptible agent
    % 2 = exposed agent
    % 3 = infected agent
    % 4 = recovered agent

    % Loop through number of susceptible agents to populate grid with
    susceptibleAgentCounter = 0;
    while (susceptibleAgentCounter < numOfSusceptible)
        % Get random coordinates
        x_coord = randi(length(cityGrid));
        y_coord = randi(length(cityGrid));
        if grid{y_coord, x_coord} == 0 && cityGrid{y_coord, x_coord} ~= -1  % Make sure not spawning on another agent or in water
            grid{y_coord, x_coord} = 1; % Spawn agent on cell
            susceptibleAgentCounter = susceptibleAgentCounter + 1;
        end
    end

    % Loop through number of infected agents to populate grid with
    infectedAgentCounter = 0;
    while (infectedAgentCounter < numOfInfected)
        % Get random coordinated
        x_coord = randi(length(grid));
        y_coord = randi(floor(length(grid) / 2)); % Only northern side for infected
        if grid{y_coord, x_coord} == 0 && cityGrid{y_coord, x_coord} ~= -1  % Make sure not spawning on another agent or in water
            grid{y_coord, x_coord} = 3; % Spawn agent on cell
            infectedAgentCounter = infectedAgentCounter + 1;
        end
    end
    gridAfterPopulation = grid; % Record and update grid
end

% Randomly moves the grid of agents
% agentGrid = referenced grid of agent cells
% cityGrid = city cells to condition if agents can move.
function agentGridAfterMove = randomlyMoveAgents(agentGrid, cityGrid)
    agentGridAfterMove = agentGrid; % Copy grid
    % Loop through each cell
    for y_coord = 1:length(agentGrid)
        for x_coord = 1:length(agentGrid)
            if agentGrid{y_coord, x_coord} > 0 % This cell is an agent
                % Pick a random direction
                randomDirection = randi(4);
                if randomDirection == 1 && y_coord ~= 1 % Going Up (1)
                    if agentGridAfterMove{y_coord - 1, x_coord} == 0 && cityGrid{y_coord - 1, x_coord} > -1 % Check if cell is valid to move to
                        % Move agent to new cell
                        agentGridAfterMove{y_coord - 1, x_coord} = agentGrid{y_coord, x_coord};
                        agentGridAfterMove{y_coord, x_coord} = 0;
                    end
                elseif randomDirection == 2 && x_coord ~= length(agentGrid) % Going Right (2)
                    if agentGridAfterMove{y_coord, x_coord + 1} == 0 && cityGrid{y_coord, x_coord + 1} > -1 % Check if cell is valid to move to
                        % Move agent to new cell
                        agentGridAfterMove{y_coord, x_coord + 1} = agentGrid{y_coord, x_coord};
                        agentGridAfterMove{y_coord, x_coord} = 0;
                    end
                elseif randomDirection == 3 && y_coord ~= length(agentGrid) % Going Down (3)
                    if agentGridAfterMove{y_coord + 1, x_coord} == 0 && cityGrid{y_coord + 1, x_coord} > -1 % Check if cell is valid to move to
                        % Move agent to new cell
                        agentGridAfterMove{y_coord + 1, x_coord} = agentGrid{y_coord, x_coord};
                        agentGridAfterMove{y_coord, x_coord} = 0;
                    end
                elseif randomDirection == 4 && x_coord ~= 1 % Going Left (4)
                    if agentGridAfterMove{y_coord, x_coord - 1} == 0 && cityGrid{y_coord, x_coord - 1} > -1 % Check if cell is valid to move to
                        % Move agent to new cell
                        agentGridAfterMove{y_coord, x_coord - 1} = agentGrid{y_coord, x_coord};
                        agentGridAfterMove{y_coord, x_coord} = 0;
                    end
                end
            end
        end
    end
end

% Goes through each agent and updates their infection state
% agentGrid = referenced grid of agent cells
% exposeChance = chance to expose an agent
% infectionChance = chance for exposed agent to be infected
% recoveryChance = chance for infected agent to recover
function agentGridAfterUpdate = updateAgentStates(agentGrid, exposeChance, infectionChance, recoveryChance)
    agentGridAfterUpdate = agentGrid; % Copy grid
    % Loop through each cell
    for y_coord = 1:length(agentGrid)
        for x_coord = 1:length(agentGrid)
            agent = agentGrid{y_coord, x_coord}; % Get agent state
            % Susceptible Agents and Recovered Agents don't do anything more.
            if agent == 2 % Exposed Agent
                agentGridAfterUpdate{y_coord, x_coord} = updateAgent(agent, infectionChance); % Attempt updating state  to infected
            elseif agent == 3 % Infected Agent
                randomDirection = randi(4); % Pick a random direction to attempt infection
                if randomDirection == 1 && y_coord ~= 1 % Picking Up (1)
                    if agentGrid{y_coord - 1, x_coord} == 1 % If there is a susceptible agent there
                         % Attempt updating state to exposed
                        exposedAgent = agentGrid{y_coord - 1, x_coord};
                        agentGridAfterUpdate{y_coord - 1, x_coord} = updateAgent(exposedAgent, exposeChance);
                    end
                elseif randomDirection == 2 && x_coord ~= length(agentGrid) % Picking Right (2)
                    if agentGrid{y_coord, x_coord + 1} == 1 % If there is a susceptible agent there
                         % Attempt updating state to exposed
                        exposedAgent = agentGrid{y_coord, x_coord + 1};
                        agentGridAfterUpdate{y_coord, x_coord + 1} = updateAgent(exposedAgent, exposeChance);
                    end
                elseif randomDirection == 3 && y_coord ~= length(agentGrid) % Picking Down (3)
                    if agentGrid{y_coord + 1, x_coord} == 1 % If there is a susceptible agent there
                         % Attempt updating state to exposed
                        exposedAgent = agentGrid{y_coord + 1, x_coord};
                        agentGridAfterUpdate{y_coord + 1, x_coord} = updateAgent(exposedAgent, exposeChance);
                    end
                elseif randomDirection == 4 && x_coord ~= 1 % Picking Left (4)
                    if agentGrid{y_coord, x_coord - 1} == 1 % If there is a susceptible agent there
                         % Attempt updating state to exposed
                        exposedAgent = agentGrid{y_coord, x_coord - 1};
                        agentGridAfterUpdate{y_coord, x_coord - 1} = updateAgent(exposedAgent, exposeChance);
                    end
                end
                agentGridAfterUpdate{y_coord, x_coord} = updateAgent(agent, recoveryChance); % Attempt to update state to recovered
            end
        end
    end
end

% Run through the chance of updating an agent to it's next state
% agent = the referenced agent
% chanceToUpdate = the chance for the agent to update to next state
function updatedAgent = updateAgent(agent, chanceToUpdate)
    updatedAgent = agent;
    if rand() <= chanceToUpdate % Run chance to update agent state
        updatedAgent = agent + 1;
    end
end

% Create a video to write to
% frameRate = framerate of the video
% index = index of the subplot for saving
% saveFileDir = the directory of the save file
function video = createVideo(frameRate, index, saveFileDir)
    vid_title = sprintf(saveFileDir + "\\video-%d.mp4", index);
    video = VideoWriter(vid_title, 'MPEG-4'); % Create a video for displaying simulation.
    video.FrameRate = frameRate;
end

% Creates a frame for the video
% cityGrid = matrix for city cells
% agentGrid = matrix for agent cells
% sizeMultiplier = increases size of frame
function colourmap = createFrame(cityGrid, agentGrid, sizeMultiplier)
    colourmap = zeros(length(cityGrid) * sizeMultiplier, length(cityGrid) * sizeMultiplier, 3); % Create empty colourmap for frame
     % Loop through each pixel
    for x_coord = 1:length(cityGrid) * sizeMultiplier
        for y_coord = 1:length(cityGrid) * sizeMultiplier
            % Calculations to make sure pixel matches with downscaled cell 
            % locations
            cityCell = cityGrid{ceil(y_coord / sizeMultiplier), ceil(x_coord / sizeMultiplier)};
            agentCell = agentGrid{ceil(y_coord / sizeMultiplier), ceil(x_coord / sizeMultiplier)};
            if agentCell == 0 % If there's no agent on cell
                if cityCell == -1 % Water Tile - Blue
                    colourmap(y_coord, x_coord, 3) = 1;
                elseif cityCell == 0 % City Tile - Grey
                    colourmap(y_coord, x_coord, 1) = 0.6;
                    colourmap(y_coord, x_coord, 2) = 0.6;
                    colourmap(y_coord, x_coord, 3) = 0.6;
                elseif cityCell == 1 % Bridge Tile - Brown
                    colourmap(y_coord, x_coord, 1) = 0.5;
                    colourmap(y_coord, x_coord, 2) = 0.3;
                    colourmap(y_coord, x_coord, 3) = 0;
                end
            else % Susceptible Agent - Black
                if agentCell == 2 % Exposed Agent - Red
                    colourmap(y_coord, x_coord, 1) = 1;
                    colourmap(y_coord, x_coord, 2) = 0;
                    colourmap(y_coord, x_coord, 3) = 0;
                elseif agentCell == 3 % Infected Agent - Green
                    colourmap(y_coord, x_coord, 1) = 0;
                    colourmap(y_coord, x_coord, 2) = 1;
                    colourmap(y_coord, x_coord, 3) = 0;
                elseif agentCell == 4 % Recovered Agent - White
                    colourmap(y_coord, x_coord, 1) = 1;
                    colourmap(y_coord, x_coord, 2) = 1;
                    colourmap(y_coord, x_coord, 3) = 1;
                end
            end
        end
    end
end

% Get number of agents for specified state
% agentGrid = grid of agents checking against
% state = specified state of agents counted
function numOfAgents = countAllAgents(agentGrid, state)
    numOfAgents = 0;
    % Loop through each cell
    for x = 1:length(agentGrid)
        for y = 1:length(agentGrid)
            if agentGrid{x, y} == state % If the agent is the state we want
                numOfAgents = numOfAgents + 1; % Add to the count
            end
        end
    end
end

% Creates the plot of each agent count
% susceptibleAgentData = tracked data of number of susceptible agents
% exposedAgentData = tracked data of number of exposed agents
% infectedAgentData = tracked data of number of infected agents
% recoveredAgentData = tracked data of number of recovered agents
% numOfBridges = number of bridges (to label each plot)
function createPlot(susceptibleAgentData, exposedAgentData, infectedAgentData, recoveredAgentData, numOfBridges)
    plotLength = length(susceptibleAgentData); % Calculate plot length
    plotHeight = 1000; % Plot Height
    x_axis = 1:plotLength; % X-Axis data
    plot(x_axis, susceptibleAgentData, 'k', ...
        x_axis, exposedAgentData, 'r', ...
        x_axis, infectedAgentData, 'g', ...
        x_axis, recoveredAgentData, 'b'); % Create plot with agent data
    % Susceptible = Black
    % Exposed = Red
    % Infected = Green
    % Recovered = Blue (Can't be white because of background)

    % Plot Settings
    xlabel('Step');
    ylabel('Num Of Agents');
    xlim([0 plotLength]);
    ylim([0 plotHeight]);
    plotTitle = sprintf('%d Bridges', numOfBridges);
    title(plotTitle);
end

% Debug function to print agents state numbers
function printAgentDistribution(agentGrid)
    fprintf("Agent Distribution: %d | %d | %d | %d\n", countAllAgents(agentGrid, 1), countAllAgents(agentGrid, 2), countAllAgents(agentGrid, 3), countAllAgents(agentGrid, 4))
end
