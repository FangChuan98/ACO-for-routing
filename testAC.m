%% ACO for routing DEMO
% Hanze 2023
% This is a simple demo for ACO-based routing. field_node.mat is a node
% data matrix storing the x-y coordinates and PV (potential value) for
% nodes. PV, to some extent, can be interpreted as the priority in routing.
clear
clc
% Initialization
numNodes = 400; % num of nodes
numAnts = 50; % num of ants
maxIterations = 100; % maximum iterations
alpha = 1; % pheromone importance coefficient
beta = 2; % heuristic information importance coefficient
rho = 0.1; % pheromone evaporation rate
Q = 100; % pheromone intensity
maxLength = 300; % maximum path length
sourceNode = 1;
targetNode = 2;

% load nodes data
load('field_node.mat');
fieldsToKeep = {'id', 'x', 'y', 'nebSet', 'pv'};
nodes = rmfield(node, setdiff(fieldnames(node), fieldsToKeep));

% Initialize the pheromone matrix
pheromoneLevels = ones(numNodes, numNodes);
bestLengths = zeros(maxIterations, 1); % Record the optimal path length for each iteration

% main loop
for iteration = 1:maxIterations
    allPaths = cell(numAnts, 2); % path store for all ants
    pathLengths = zeros(numAnts, 1); % storage path length

    for ant = 1:numAnts
        tabuList = zeros(1, numNodes); % tabu list
        currentNode = sourceNode; 
        path = currentNode; 
        tabuList(currentNode) = 1; 
        success = false; 

        while ~isempty(currentNode) && currentNode ~= targetNode && length(path) < maxLength
            if ismember(targetNode, nodes(currentNode).nebSet)
                nextNode = targetNode;
                path = [path, nextNode];
                currentNode = nextNode;
                success = true;
            else
                neighbors = nodes(currentNode).nebSet;
                neighbors = neighbors(~ismember(neighbors, find(tabuList)));
                probabilities = zeros(1, numel(neighbors));
                for i = 1:numel(neighbors)
                    nextNode = neighbors(i);
                    tau = pheromoneLevels(currentNode, nextNode)^alpha;
                    eta_dist = 1 / distance(nodes(nextNode), nodes(targetNode));
                    eta_pv = (1 / nodes(nextNode).pv)^beta;
                    probabilities(i) = tau * eta_dist * eta_pv;
                end
                probabilities = probabilities / sum(probabilities);
                nextNodeIndex = find(rand <= cumsum(probabilities), 1);
                nextNode = neighbors(nextNodeIndex);

                path = [path, nextNode];
                tabuList(nextNode) = 1;
                currentNode = nextNode;
            end
        end

        % save to all Paths
        allPaths{ant, 1} = path;
        allPaths{ant, 2} = success;
        if success
            pathLengths(ant) = calculatePathLength(path, nodes);
        else
            pathLengths(ant) = Inf;
        end
    end

    % updata
    deltaPheromone = zeros(numNodes, numNodes);
    for ant = 1:numAnts
        path = allPaths{ant, 1};
        if allPaths{ant, 2} % if reach
            L = pathLengths(ant);
            for i = 1:(length(path)-1)
                deltaPheromone(path(i), path(i+1)) = deltaPheromone(path(i), path(i+1)) + Q / L;
            end
        end
    end
    pheromoneLevels = (1 - rho) * pheromoneLevels + deltaPheromone;

    % store best path lengths
    [bestLengths(iteration), idx] = min(pathLengths);
    bestPaths{iteration} = allPaths{idx, 1};
end

% draw curve
figure,
plot(bestLengths);
xlabel('Iteration');
ylabel('Best Path Length');
title('Convergence of ACO Algorithm');

% draw best path
[~, bestIdx] = min(bestLengths);
bestPath = bestPaths{bestIdx};
figure;
hold on;
for i = 1:length(bestPath)-1
    node1 = nodes(bestPath(i));
    node2 = nodes(bestPath(i+1));
    plot([node1.x, node2.x], [node1.y, node2.y], 'r-o'); 
end

for i = 1:numNodes
    plot(nodes(i).x, nodes(i).y, 'b-o'); 
end

text(nodes(sourceNode).x, nodes(sourceNode).y, 'Source');
text(nodes(targetNode).x, nodes(targetNode).y, 'Target');

xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Best Path Found by ACO');
hold off;

function len = calculatePathLength(path, nodes)
    len = 0;
    for i = 1:(length(path)-1)
        len = len + distance(nodes(path(i)), nodes(path(i+1)));
    end
end

function dist = distance(node1, node2)
    dist = sqrt((node1.x - node2.x)^2 + (node1.y - node2.y)^2);
end
