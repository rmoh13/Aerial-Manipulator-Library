%{
inputs are the waypoints and associated initial conditions/values, and
time segments. 

output is position, velocity, acceleration, jerk at any point in time

assume that we're working with minimum snap trajectory, and that
we're given the starting and end position, velocity, acceleration, and jerk
if we weren't given the starting and end values, that would make things 
much more difficult
so this means assume we're working with polynomial trajectories that are
7th order polynomials
%}

%{ 
waypoints is a cell array matrix, 
each cell is a column vector, and 
the first cell is assumed to be a column vector of all the positions,
the next columns are column vectors of vel, accel, jerk at each position
in that order

times is a column vector of time points for the time segments
assume they're given in order in the waypoints cell and if there's 
nothing at a cer
%}

function vals = multiSegmentPlanner(waypoints, times)
    % in Ax = b, vals is x which is the constraints column vector
    % make waypoints rows, assume each position is given
    numT = size(times, 1); %#ok<NASGU>
    positions = waypoints{1,1};
    numWP = size(positions, 1);
    % this should be true: numWP == numT
    % A = zeros(8 * (numWP - 1), 8 * (numWP - 1));
    A = zeros(2 * (numWP - 2), 8 * (numWP - 1));
    b = zeros(2 * numWP - 2, 1);

    
    % fill up b w/position constraints
    b(1,1) = positions(1,1);
    b(size(b, 1), 1) = positions(numWP, 1);
    for j = 2:size(b,1)-1
        if mod(j, 2) == 0
            i = (j / 2) + 1;
            b(j, 1) = positions(i, 1);
            b(j + 1, 1) = positions(i, 1);
        end
    end
    
    % fill up A w/position constraints
    % get the first and last position constraint parts of A
    A(1, 1 : 8) = ...
            [times(1)^7 times(1)^6 times(1)^5 times(1)^4 times(1)^3 ...
            times(1)^2 times(1)^1 times(1)^0];
        
    A(size(b,1), ...
    8*(numWP - 1) - 7 : 8*(numWP - 1)) = ...
            [times(numWP,1)^7 times(numWP,1)^6 times(numWP,1)^5 ...
            times(numWP,1)^4 times(numWP,1)^3 ...
            times(numWP,1)^2 times(numWP,1)^1 times(numWP,1)^0];
    
    % iterate in between
    for n = 2:size(b,1)-2
        if mod(n, 2) == 0
            i = (n / 2) + 1;
            A(n, (8*ceil((n-2) / 2)) + 1 : (8*ceil((n-2) / 2) + 8)) = ...
             [times(i,1)^7 times(i,1)^6 times(i,1)^5 times(i,1)^4 times(i,1)^3 ...
             times(i,1)^2 times(i,1)^1 times(i,1)^0];
            A(n + 1, (8*ceil((n-1) / 2)) + 1 : (8*ceil((n-1) / 2) + 8)) = ...
             [times(i,1)^7 times(i,1)^6 times(i,1)^5 times(i,1)^4 ...
             times(i,1)^3 times(i,1)^2 times(i,1)^1 times(i,1)^0];
        end
    end
    
    %{ 
    iterate over the rest of the waypoints cell array, each of which is
    a column vector of the velocities, accelerations, and jerks at
    each waypoint
    %}
    
    numTraj = numWP - 1; %#ok<NASGU>
    % starting row to add the value constraints to A
    % addValsA = size(b,1) + 1;
    
    % bVals is the column vec of value constraints we will combine b with
    % can't really precompute here, so start w/ first val, iterate to last
    bVals = [waypoints{1,2}];
    AValsFirstRow = [];
    % do the first row of the corresponding value constraints in A
    for i = 1:size(waypoints{1, 2}, 1)
        if i == 1
            AValsFirstRow = [AValsFirstRow ; calcDerTraj(1, 1) ,  zeros(1, size(A, 2) - size(calcDerTraj(1, 1), 2))]; 
        elseif i == 2
            AValsFirstRow = [AValsFirstRow ; calcDerTraj(1, 2) ,  zeros(1, size(A, 2) - size(calcDerTraj(1, 2), 2))]; 
        elseif i == 3
            AValsFirstRow = [AValsFirstRow ; calcDerTraj(1, 3) ,  zeros(1, size(A, 2) - size(calcDerTraj(1, 3), 2))]; 
        end
    end
    A = [A ; AValsFirstRow];
        
    % now go in between
    j = 3;
    while j < size(waypoints, 2)
        bVals = [bVals ; waypoints{1,j}; waypoints{1,j}]; 
        j = j + 1;
    end
    bVals = [bVals ; waypoints{1, size(waypoints, 2)}];
    b = [b ; bVals]; %#ok<NASGU>
    AInBetween = [];
    for i = 3:size(waypoints, 2) - 1
        timeIndex = i - 1;
        for k = 1:size(waypoints{1, i}, 1)
             if k == 1
                 der1Traj = calcDerTraj(timeIndex, 1);
                 currRow = zeros(1, 8 * (numWP - 1))
                 currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der1Traj, 2)) ...
                 = der1Traj
                AInBetween = [AInBetween ; currRow] %#ok<*AGROW>
             elseif k == 2
                 der2Traj = calcDerTraj(timeIndex, 2);
                 currRow = zeros(1, 8 * (numWP - 1))
                 currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der2Traj, 2)) ...
                 = der2Traj
                AInBetween = [AInBetween ; currRow] %#ok<*AGROW>
             elseif k == 3
                 der3Traj = calcDerTraj(timeIndex, 3);
                 currRow = zeros(1, 8 * (numWP - 1))
                 currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der3Traj, 2)) ...
                 = der3Traj
                AInBetween = [AInBetween ; currRow] %#ok<*AGROW>
             end         
        end
        for k = 1:size(waypoints{1, i}, 1)
             if k == 1
                 der1Traj = calcDerTraj(timeIndex, 1);
                 currRow = zeros(1, 8 * (numWP - 1))
                 currRow(1, 8*(i-2) + 1 : 8*(i-2) + size(der1Traj, 2)) ...
                 = der1Traj
                AInBetween = [AInBetween ; currRow] %#ok<*AGROW>
             elseif k == 2
                 der2Traj = calcDerTraj(timeIndex, 2);
                 currRow = zeros(1, 8 * (numWP - 1))
                 currRow(1, 8*(i-2) + 1 : 8*(i-2) + size(der2Traj, 2)) ...
                 = der2Traj
                AInBetween = [AInBetween ; currRow] %#ok<*AGROW>
             elseif k == 3
                 der3Traj = calcDerTraj(timeIndex, 3);
                 currRow = zeros(1, 8 * (numWP - 1))
                 currRow(1, 8*(i-2) + 1 : 8*(i-2) + size(der3Traj, 2)) ...
                 = der3Traj
                AInBetween = [AInBetween ; currRow] %#ok<*AGROW>
             end         
        end
    end
    A = [A ; AInBetween];
    % now do the last row
    AValsLastRow = [];
    for i = 1:size(waypoints{1, size(waypoints, 2)}, 1)
        if i == 1
            AValsLastRow = ...
            [AValsLastRow ; zeros(1, size(A, 2) - size(calcDerTraj(size(waypoints, 2) - 1, 1), 2)), ...
            calcDerTraj(size(waypoints, 2) - 1, 1)]; 
        elseif i == 2
            AValsLastRow = ...
            [AValsLastRow ; zeros(1, size(A, 2) - size(calcDerTraj(size(waypoints, 2) - 1, 2), 2)), ...
            calcDerTraj(size(waypoints, 2) - 1, 2)]; 
        elseif i == 3
            AValsLastRow = ...
            [AValsLastRow ; zeros(1, size(A, 2) - size(calcDerTraj(size(waypoints, 2) - 1, 3), 2)), ...
            calcDerTraj(size(waypoints, 2) - 1, 3)]; 
        end
    end
    A = [A ; AValsLastRow];
    
    vals = A;
end

function derTraj = calcDerTraj(timeIndex, nthDer)
    % we are assuming a 7th order polynomial for the trajectory (min. snap)
    traj = [times(timeIndex, 1)^7 times(timeIndex, 1)^6 ...
            times(timeIndex, 1)^5 times(timeIndex, 1)^4 ...
            times(timeIndex, 1)^3 times(timeIndex, 1)^2 ...
            times(timeIndex, 1)^1 times(timeIndex, 1)^0];
        sizeTraj = size(traj,2);
        dTraj = (sizeTraj-1:-1:1).*traj(1:sizeTraj-1);
        size2Traj = size(dTraj,2);
        d2Traj = (size2Traj-1:-1:1).*dTraj(1:size2Traj-1);
        size3Traj = size(d2Traj,2);
        d3Traj = (size3Traj-1:-1:1).*d2Traj(1:size3Traj-1);
        if nthDer == 0
            derTraj = traj;
        elseif nthDer == 1
            derTraj = dTraj;
        elseif nthDer == 2
            derTraj = d2Traj;
        elseif nthDer == 3
            derTraj = d3Traj;
        end
end