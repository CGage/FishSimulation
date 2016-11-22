function [] = itsAFishEatFishWorld()
    close all;     
    width = 20;
    height = 20;
    timeSteps = 1000;
    currentSpeed = 1;
    maxFishPerCell = 5; % This value must also be changed in the trim edges function.
    maxPlanktonPerStep = 100;
    initPlankton =maxPlanktonPerStep * height;
    
    deltax = 1; deltay = 1;
    
    CA = randi(maxFishPerCell, height, width);
    RW = randi(height, 2, initPlankton);
    cmap = makeCmap(maxFishPerCell);
    for i = 1: timeSteps
         currentSpeed = randi(2); % Simulates a changing current speed
        planktonCount1 = countPlankton(RW, height, width);
        plotFishies(RW, CA, width, height, cmap, planktonCount1);
        pause(0.5);     
        RW = addPlankton(RW, height, maxPlanktonPerStep);
    	RW = computeRW(RW, currentSpeed, width, height, deltax, deltay);
        planktonCount2 = countPlankton(RW, height, width);
%         CA = migrate(planktonCount1, planktonCount2, CA, height, width);
        [CA, RW] = eatPlankton(CA, RW, planktonCount2, maxFishPerCell, width, height);        
        CA = computeCA( CA, maxFishPerCell, width, height);        
    end
end


% Adds a random number of plankton to the random walk simulation
function RW = addPlankton(RW, height, maxPlanktonPerStep)
    newPlankton = ones(2, randi(maxPlanktonPerStep));
    newPlankton(2, :) = randi([0 height], [1 size(newPlankton,2)]);    
    RW = [RW newPlankton];
end

% Computes the random walk simulation for one time step and removes any
% plankton that leave the area.
function RW = computeRW(RW, currentSpeed, width, height, deltax, deltay)
    RW(1, :) = RW(1, :) + currentSpeed; % simulates the flow of the water

    r = rand(1,size(RW, 2)); % generate M random numbers between 0 and 1
    
    left_mask  = r < 0.25;                            % mask identifying the left-moving particles
    RW(1, left_mask) = RW(1, left_mask)  - deltax;    % move those particles left
    RW(2, left_mask) = RW(2, left_mask);

    right_mask = r >= 0.25 & r < 0.5;                 % mask identifying the right-moving particles
    RW(1, right_mask) = RW(1, right_mask) + deltax;   % move those particles right
    RW(2, right_mask) = RW(2, right_mask);

    down_mask = r >= 0.5 & r < 0.75;                  % mask identifying the down-moving particles
    RW(1, down_mask) = RW(1, down_mask);   
    RW(2, down_mask) = RW(2, down_mask) - deltay;     % move those particles down

    up_mask = r >= 0.75;                              % mask identifying the up-moving particles
    RW(1, up_mask) = RW(1, up_mask);   
    RW(2, up_mask) = RW(2, up_mask) + deltay;         % move those particles up

    outOfRangeMask = RW(1, :) > width | RW(2, :) >= height | RW(2, :) <= 0;
    RW(:, outOfRangeMask) = [];
end 

% Returns a 2d array containing the number of plankton present in each cell
function planktonCount = countPlankton(RW, height, width)
    planktonCount = zeros([height width]);
    for i = 1: size(RW, 2)
        planktonCount(RW(2, i), RW(1, i)) = planktonCount(RW(2, i), RW(1, i)) + 1;
    end
end


% Enforces the minimum amount of fish allowed in a cell.
function cell =trimEdges(cell)
    if (cell < 0)
        cell = 0;
    end
end
% Ensures the fish dont consume more food than possible
function cell = foodEaten(pk, mep)
    if (pk > mep)
        cell = mep;
    else
       cell = pk; 
    end       
end

function [CA, RW] = eatPlankton(CA, RW, planktonCount, maxFishPerCell, width, height)   
    lives = (CA >= 2 & planktonCount > 0);
    dies = (planktonCount < 1);
    mep = lives .* (maxFishPerCell - CA); % maximum plankton that can be eaten.
    planktonEaten = arrayfun(@(x, y) foodEaten(x, y),planktonCount , mep);
    CA = CA + planktonEaten - dies;
    CA = arrayfun(@(x)trimEdges(x), CA);
    planktonEaten = planktonEaten + (CA == maxFishPerCell & planktonCount > 0); % If there are maxFishPerCell fish in the cell they consume one plankton and stay at the same number.
    for h = 1: width
        for j = 1: height
            for k = 1: planktonEaten(h, j)
                for l = size(RW, 2): -1: 1
                    if(RW(1, l) == j && RW(2, l) == h)
                        RW(:, l) = []; % Removes a plankton from the random walk when it is consumed.
                        break;
                    end
                    
                end
            end
        end
    end
end

function [Max, index] = elementwiseMax(n, s, e, w, ne, nw, se, sw, c)
[Max index] = max([n, s, e, w, ne, nw, se, sw, c]);

end
function CA = migrate(planktonCount1, planktonCount2, CA, height, width)
    planktonDiff = planktonCount2 - planktonCount1;
    north = [height 1:height-1]; % indices of north neighbour
    east  = [2:width 1];       % indices of east neighbour
    south = [2:height 1];       % indices of south neighbour
    west  = [width 1:width-1];     % indices of west neighbour
%     [bestMove index] = max([planktonDiff(north, :) planktonDiff(south, :) planktonDiff(:, east)...
%         planktonDiff(:, west) planktonDiff(north, east) planktonDiff(north, west)...
%         planktonDiff(south, east) planktonDiff(south, west)])
    n = planktonDiff(north, :);
    s = planktonDiff(south, :);
    e = planktonDiff(:, east);
    w = planktonDiff(:, west);
    ne = planktonDiff(north, east);
    nw = planktonDiff(north, west);
    se = planktonDiff(south, east);
    sw = planktonDiff(south, west);
    c = planktonDiff();
    [best, index] = arrayfun(@(x1, x2, x3, x4, x5, x6, x7, x8, x9) ...
       elementwiseMax(x1, x2, x3, x4, x5, x6, x7, x8, x9) ...
       ,n, s, e, w, ne, nw, se, sw, c);

    pause(10)
    
end

function [CA] = computeCA(CA, maxFishPerCell, width, height) 
    north = [height 1:height-1]; % indices of north neighbour
    east  = [2:width 1];       % indices of east neighbour
    south = [2:height 1];       % indices of south neighbour
    west  = [width 1:width-1];     % indices of west neighbour


    % Count how many live neighbours each cell has in its Moore neighbourhood
    
    live_neighbours = CA(north, :) + CA(south, :) + CA(:, east) + CA(:, west) ...
    	+ CA(north, east) + CA(north, west) + CA(south, east) + CA(south, west);              
    % There are only 2 ways that a cell can live in the Game of Life:
    alive_rule_1 = live_neighbours <= (maxFishPerCell*4) & CA ~= 0;        % a cell lives if it has 3 live neighbours
    alive_rule_2 = live_neighbours >= 2 & live_neighbours <= (maxFishPerCell * 4) & CA == 0; % a cell lives if it's alive already, and has 2 live neighbours
    fishBornInNewCell = 2;
    % These two rules determine the new state of every element
    CA = (CA .* alive_rule_1) + (fishBornInNewCell .* alive_rule_2);
    % 
    
end

% Returns a colourmap scaled over the maximum number of possible fish 
% allowed in each cell.
function cmap = makeCmap(maxFishPerCell)   
    cmap = zeros([maxFishPerCell 3]);
    row1 = linspace(0.58, 1, maxFishPerCell)';
    row2 = linspace(0.83, 1, maxFishPerCell)';
    row3 = linspace(0.83, 1, maxFishPerCell)';
    cmap(:, 1) = row1;
    cmap(:, 2) = row2;
    cmap(:, 3) = row3;
    cmap = flipud(cmap);
end

% Plots both the cellular automata and the random walk on the same figure.
function [] = plotFishies(RW, CA, width, height, cmap, planktonCount)
    hold off;
    imshow(CA, cmap,'InitialMagnification','fit');
    hold on;

    % This plot is suitable for larger problem sizes.
    scatter(RW(1, :), RW(2, :), '.k'); 

    % Warning: Plotting numbers onto the figure is computationally expensive so
    % it is wise to keep the simulations sizes small when using this 
    % visualisation "mode", ie under 30 width *30 height.
%     y = linspace(1, height, height); 
%     for i = 1: height                
%         x = linspace(i, i, width);  
%         text(x, y ,num2str(planktonCount(:, i))); 
%     end
end