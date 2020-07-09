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
%}

% waypoints must be a cell array, times must be a column vector
%{
note that for any value in the intermediary WPs, if the user 
can't/doesn't want to give a value at a specific WP, either leave
the entire column vector empty so they don't provide any ICs at the WP
or put a pi for DNE values. Example, if user wants to only provide
acceleration at one of the waypoints, give a value for acceleration
at that space in the column and put pi everywhere else
%}
function vals = multiSegmentPlanner(waypoints, times)
    % in Ax = b, vals is x which is the constraints column vector
    % make waypoints rows, assume each position is given
    numT = size(times, 1);
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
            AValsFirstRow = [AValsFirstRow ; calcDerTraj(times, 1, 1) ,  zeros(1, size(A, 2) - size(calcDerTraj(times, 1, 1), 2))]; 
        elseif i == 2
            AValsFirstRow = [AValsFirstRow ; calcDerTraj(times, 1, 2) ,  zeros(1, size(A, 2) - size(calcDerTraj(times, 1, 2), 2))]; 
        elseif i == 3
            AValsFirstRow = [AValsFirstRow ; calcDerTraj(times, 1, 3) ,  zeros(1, size(A, 2) - size(calcDerTraj(times, 1, 3), 2))]; 
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
    b = [b ; bVals];
    AInBetween = [];
    for i = 3:size(waypoints, 2) - 1
        timeIndex = i - 1;
        for k = 1:size(waypoints{1, i}, 1)
             if k == 1
                 der1Traj = calcDerTraj(times, timeIndex, 1);
                 currRow = zeros(1, 8 * (numWP - 1));
                 currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der1Traj, 2)) ...
                 = der1Traj;
                 AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
             elseif k == 2
                 der2Traj = calcDerTraj(times, timeIndex, 2);
                 currRow = zeros(1, 8 * (numWP - 1));
                 currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der2Traj, 2)) ...
                 = der2Traj;
                 AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
             elseif k == 3
                 der3Traj = calcDerTraj(times, timeIndex, 3);
                 currRow = zeros(1, 8 * (numWP - 1));
                 currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der3Traj, 2)) ...
                 = der3Traj;
                 AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
             end         
        end
        for k = 1:size(waypoints{1, i}, 1)
             if k == 1
                 der1Traj = calcDerTraj(times, timeIndex, 1);
                 currRow = zeros(1, 8 * (numWP - 1));
                 currRow(1, 8*(i-2) + 1 : 8*(i-2) + size(der1Traj, 2)) ...
                 = der1Traj;
                 AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
             elseif k == 2
                 der2Traj = calcDerTraj(times, timeIndex, 2);
                 currRow = zeros(1, 8 * (numWP - 1));
                 currRow(1, 8*(i-2) + 1 : 8*(i-2) + size(der2Traj, 2)) ...
                 = der2Traj;
                 AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
             elseif k == 3
                 der3Traj = calcDerTraj(times, timeIndex, 3);
                 currRow = zeros(1, 8 * (numWP - 1));
                 currRow(1, 8*(i-2) + 1 : 8*(i-2) + size(der3Traj, 2)) ...
                 = der3Traj;
                 AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
             end         
        end
    end
    A = [A ; AInBetween];
    % now do the last row
    AValsLastRow = [];
    lastWPIndex = size(waypoints, 2) - 1;
    for k = 1:size(waypoints{1, size(waypoints, 2)}, 1)
         if k == 1
             der1Traj = calcDerTraj(times, lastWPIndex, 1);
             currRow = zeros(1, 8 * (numWP - 1));
             currRow(1, 8*(lastWPIndex-2) + 1 : 8*(lastWPIndex-2) + size(der1Traj, 2)) ...
             = der1Traj;
             AValsLastRow = [AValsLastRow ; currRow]; %#ok<*AGROW>
         elseif k == 2
             der2Traj = calcDerTraj(times, lastWPIndex, 2);
             currRow = zeros(1, 8 * (numWP - 1));
             currRow(1, 8*(lastWPIndex-2) + 1 : 8*(lastWPIndex-2) + size(der2Traj, 2)) ...
             = der2Traj;
             AValsLastRow = [AValsLastRow ; currRow]; %#ok<*AGROW>
         elseif k == 3
             der3Traj = calcDerTraj(times, lastWPIndex, 3);
             currRow = zeros(1, 8 * (numWP - 1));
             currRow(1, 8*(lastWPIndex-2) + 1 : 8*(lastWPIndex-2) + size(der3Traj, 2)) ...
             = der3Traj;
             AValsLastRow = [AValsLastRow ; currRow]; %#ok<*AGROW>
         end         
    end
    A = [A ; AValsLastRow];
    
    %{ 
    the above is working for when all values, meaning the velocities,
    accelerations, and jerks are given at each WP. Now, we have to 
    account for values not given and go case by case to give the
    linear equality constraints as such in the matrix A
    %}
    
    % assume ICs at starting and ending WPs must be given
    for i = 3:size(waypoints, 2) - 1
        timeIndex = i - 1;
        indicesUnknowns = find(waypoints{1, i} == pi);
        indicesKnowns = find(waypoints{1, i} ~= pi); %#ok<NASGU>
        for k = 1:size(indicesUnknowns, 2)
            currentIndexUnknown = indicesUnknown(1, k);
            waypoints{1,i}(currentIndexUnknown, 1);
             if k == 1
                 der1Traj = calcDerTraj(times, timeIndex, 1);
                 currRow = zeros(1, 8 * (numWP - 1));
                 currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der1Traj, 2)) ...
                 = der1Traj;
                 AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
             elseif k == 2
                 der2Traj = calcDerTraj(times, timeIndex, 2);
                 currRow = zeros(1, 8 * (numWP - 1));
                 currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der2Traj, 2)) ...
                 = der2Traj;
                 AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
             elseif k == 3
                 der3Traj = calcDerTraj(times, timeIndex, 3);
                 currRow = zeros(1, 8 * (numWP - 1));
                 currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der3Traj, 2)) ...
                 = der3Traj;
                 AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
             end         
        end
    end
    
    x = A \ b;
    sizeSol = size(x, 1);
    
    
    % now we do the plotting
    i = 1;
    timeIndex = 1;
    subplot(4,1,1);
    hold on
    while i < sizeSol
        posTraj = x(i:i+7).';
        t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
        y = posTraj;
        plot(t, polyval(y,t))
        title("position vs time")
        xlabel("t")
        ylabel("x(t)")
        i = i + 8;
        timeIndex = timeIndex + 1;
    end
    hold off
    
    i = 1;
    timeIndex = 1;
    subplot(4,1,2);
    hold on
    while i < sizeSol
        posTraj = x(i:i+7).';
        velTraj = polyder(posTraj);
        t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
        y = velTraj;
        plot(t, polyval(y,t))
        title("velocity vs time")
        xlabel("t")
        ylabel("x'(t)")
        i = i + 8;
        timeIndex = timeIndex + 1;
    end
    hold off
    
    i = 1;
    timeIndex = 1;
    subplot(4,1,3);
    hold on
    while i < sizeSol
        posTraj = x(i:i+7).';
        velTraj = polyder(posTraj);
        accelTraj = polyder(velTraj);
        t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
        y = accelTraj;
        plot(t, polyval(y,t))
        title("acceleration vs time")
        xlabel("t")
        ylabel("x''(t)")
        i = i + 8;
        timeIndex = timeIndex + 1;
    end
    hold off
    
    i = 1;
    timeIndex = 1;
    subplot(4,1,4);
    hold on
    while i < sizeSol
        posTraj = x(i:i+7).';
        velTraj = polyder(posTraj);
        accelTraj = polyder(velTraj);
        jerkTraj = polyder(accelTraj);
        t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
        y = jerkTraj;
        plot(t, polyval(y,t))
        title("jerk vs time")
        xlabel("t")
        ylabel("x'''(t)")
        i = i + 8;
        timeIndex = timeIndex + 1;
    end
    hold off
    
    
    vals = {A, b, x};
end

function derTraj = calcDerTraj(times, timeIndex, nthDer)
    % we are assuming a 7th order polynomial for the trajectory (min. snap)
    traj = [times(timeIndex, 1)^7 times(timeIndex, 1)^6 ...
            times(timeIndex, 1)^5 times(timeIndex, 1)^4 ...
            times(timeIndex, 1)^3 times(timeIndex, 1)^2 ...
            times(timeIndex, 1)^1 times(timeIndex, 1)^0];    
        normCoeff = [1 1 1 1 1 1 1 1];
        derCoeff = polyder(normCoeff);
        der2Coeff = polyder(derCoeff);   
        der3Coeff = polyder(der2Coeff);   
        der4Coeff = polyder(der3Coeff);
        der5Coeff = polyder(der4Coeff);
        der6Coeff = polyder(der5Coeff);
        if nthDer == 0
            derTraj = traj;
        elseif nthDer == 1
            derTraj = derCoeff.*traj(nthDer+1:size(traj, 2));
        elseif nthDer == 2
            derTraj = der2Coeff.*traj(nthDer+1:size(traj, 2));
        elseif nthDer == 3
            derTraj = der3Coeff.*traj(nthDer+1:size(traj, 2));
        elseif nthDer == 4
            derTraj = der4Coeff.*traj(nthDer+1:size(traj, 2));
        elseif nthDer == 5
            derTraj = der5Coeff.*traj(nthDer+1:size(traj, 2));
        elseif nthDer == 6
            derTraj = der6Coeff.*traj(nthDer+1:size(traj, 2));
        end
end