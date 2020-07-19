classdef MultiSegmentTrajPlanner < Trajectory
    %MULTISEGMENTTRAJPLANNER Multi-Segment Trajectory Planner
    %{  
        assume that we're working with minimum snap trajectory, and that
        we're given the starting and end position, velocity, acceleration, and jerk
        if we weren't given the starting and end values, that would make things 
        much more difficult
        so this means assume we're working with polynomial trajectories that are
        7th order polynomials
    %}
    
    properties (Access = private)
        waypoints
        times
    end
    
    methods
        function obj = MultiSegmentTrajPlanner(waypoints, times, dimensions)
            %MULTISEGMENTTRAJPLANNER Constructor
            %   inputs are the waypoints and associated initial conditions/values, and time segments. 
            %{ 
            waypoints is a cell array matrix, 
            each cell is a column vector, and 
            the first cell is assumed to be a column vector of all the positions,
            the next columns are column vectors of vel, accel, jerk at each position
            in that order

            times is a column vector of time points for the time segments
            %}

            %{
            note that for any value in the intermediary WPs, if the user 
            can't/doesn't want to give a value at a specific WP, either leave
            the entire column vector empty so they don't provide any ICs at the WP
            or put a NaN for DNE values. Example, if user wants to only provide
            acceleration at one of the waypoints, give a value for acceleration
            at that space in the column and put NaN everywhere else
            %}
            obj.waypoints = waypoints;
            obj.times = times;
            obj.duration = times(size(times, 1), 1) - times(1, 1);
            obj.dimensions = dimensions;
        end
               
        function times = get.times(obj)
            times = obj.times;
        end
        
        function waypoints = get.waypoints(obj)
            waypoints = obj.waypoints;
        end
        
        function obj = set.times(obj, newTimes)
            obj.times = newTimes;
        end
        
        function obj = set.waypoints(obj, newWaypoints)
            obj.waypoints = newWaypoints;
        end
        
        function vals = generateTrajAndCoeff(obj)
            % in Ax = b, vals is x which is the constraints column vector
            % make waypoints rows, assume each position is given
            waypoints = obj.waypoints; %#ok<*PROP>
            times = obj.times;
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
    
            % starting row to add the value constraints to A
            % addValsA = size(b,1) + 1;
    
            % bVals is the column vec of value constraints we will combine b with
            % can't really precompute here, so start w/ first val, iterate to last
            bVals = [waypoints{1,2}];
            AValsFirstRow = [];
            % do the first row of the corresponding value constraints in A
            for i = 1:size(waypoints{1, 2}, 1)
                if i == 1 && ~isnan(waypoints{1,2}(i, 1))
                    AValsFirstRow = [AValsFirstRow ; calcDerTraj(times, 1, 1) , ...  
                    zeros(1, size(A, 2) - size(calcDerTraj(times, 1, 1), 2))]; 
                elseif i == 2 && ~isnan(waypoints{1,2}(i, 1))
                    AValsFirstRow = [AValsFirstRow ; calcDerTraj(times, 1, 2) , ...
                    zeros(1, size(A, 2) - size(calcDerTraj(times, 1, 2), 2))]; 
                elseif i == 3 && ~isnan(waypoints{1,2}(i, 1))
                    AValsFirstRow = [AValsFirstRow ; calcDerTraj(times, 1, 3) , ... 
                    zeros(1, size(A, 2) - size(calcDerTraj(times, 1, 3), 2))]; 
                end
            end
            A = [A ; AValsFirstRow];
        
            % now go in between
            j = 3;
            while j < size(waypoints, 2)
                currentCol = [];
                for k = 1:size(waypoints{1,j}, 1)
                    if ~isnan(waypoints{1,j}(k, 1))
                        currentCol = [currentCol ; waypoints{1,j}(k, 1)];
                    end
                end
                bVals = [bVals ; currentCol; currentCol]; 
                j = j + 1;
            end
            bVals = [bVals ; waypoints{1, size(waypoints, 2)}];
            b = [b ; bVals];
            AInBetween = [];
            for i = 3:size(waypoints, 2) - 1
                timeIndex = i - 1;
                for k = 1:size(waypoints{1, i}, 1)
                     if k == 1 && ~isnan(waypoints{1,i}(k, 1))
                         der1Traj = calcDerTraj(times, timeIndex, 1);
                         currRow = zeros(1, 8 * (numWP - 1));
                         currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der1Traj, 2)) ...
                         = der1Traj;
                         AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
                     elseif k == 2 && ~isnan(waypoints{1,i}(k, 1))
                         der2Traj = calcDerTraj(times, timeIndex, 2);
                         currRow = zeros(1, 8 * (numWP - 1));
                         currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der2Traj, 2)) ...
                         = der2Traj;
                         AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
                     elseif k == 3 && ~isnan(waypoints{1,i}(k, 1))
                         der3Traj = calcDerTraj(times, timeIndex, 3);
                         currRow = zeros(1, 8 * (numWP - 1));
                         currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der3Traj, 2)) ...
                         = der3Traj;
                         AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
                     end         
                end
                for k = 1:size(waypoints{1, i}, 1)
                     if k == 1 && ~isnan(waypoints{1,i}(k, 1))
                         der1Traj = calcDerTraj(times, timeIndex, 1);
                         currRow = zeros(1, 8 * (numWP - 1));
                         currRow(1, 8*(i-2) + 1 : 8*(i-2) + size(der1Traj, 2)) ...
                         = der1Traj;
                         AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
                     elseif k == 2 && ~isnan(waypoints{1,i}(k, 1))
                         der2Traj = calcDerTraj(times, timeIndex, 2);
                         currRow = zeros(1, 8 * (numWP - 1));
                         currRow(1, 8*(i-2) + 1 : 8*(i-2) + size(der2Traj, 2)) ...
                         = der2Traj;
                         AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
                     elseif k == 3 && ~isnan(waypoints{1,i}(k, 1))
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
            lastWNaNndex = size(waypoints, 2) - 1;
            for k = 1:size(waypoints{1, size(waypoints, 2)}, 1)
                 if k == 1 && ~isnan(waypoints{1,size(waypoints, 2)}(k, 1))
                     der1Traj = calcDerTraj(times, lastWNaNndex, 1);
                     currRow = zeros(1, 8 * (numWP - 1));
                     currRow(1, 8*(lastWNaNndex-2) + 1 : 8*(lastWNaNndex-2) + size(der1Traj, 2)) ...
                     = der1Traj;
                     AValsLastRow = [AValsLastRow ; currRow]; %#ok<*AGROW>
                 elseif k == 2 && ~isnan(waypoints{1,size(waypoints, 2)}(k, 1))
                     der2Traj = calcDerTraj(times, lastWNaNndex, 2);
                     currRow = zeros(1, 8 * (numWP - 1));
                     currRow(1, 8*(lastWNaNndex-2) + 1 : 8*(lastWNaNndex-2) + size(der2Traj, 2)) ...
                     = der2Traj;
                     AValsLastRow = [AValsLastRow ; currRow]; %#ok<*AGROW>
                 elseif k == 3 && ~isnan(waypoints{1,size(waypoints, 2)}(k, 1))
                     der3Traj = calcDerTraj(times, lastWNaNndex, 3);
                     currRow = zeros(1, 8 * (numWP - 1));
                     currRow(1, 8*(lastWNaNndex-2) + 1 : 8*(lastWNaNndex-2) + size(der3Traj, 2)) ...
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
            AUnknowns = [];
            bUnknowns = [];
            for i = 3:size(waypoints, 2) - 1
                timeIndex = i - 1;
                indicesUnknowns = find(isnan(waypoints{1, i}));
                if size(indicesUnknowns, 1) ~= 0
                    for k = 1:size(waypoints{1, i}, 1)
                        %currentIndexUnknown = indicesUnknowns(1, k);
                        %waypoints{1,i}(currentIndexUnknown, 1);
                         %{
                         if velocity is unknown, add two constraints: the velocity
                         equality constraint and the 4th derivative constraint
                         %}
                         if k == 1 && isnan(waypoints{1,i}(k, 1))
                             % velocity equality constraint
                             der1Traj = calcDerTraj(times, timeIndex, 1);
                             currRow = zeros(1, 8 * (numWP - 1));
                             currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der1Traj, 2)) ...
                             = der1Traj;
                             currRow2 = zeros(1, 8 * (numWP - 1));
                             der1Traj2 = -1 * der1Traj;
                             currRow2(1, 8*(i-2) + 1 : 8*(i-2) + size(der1Traj, 2)) ...
                             = der1Traj2;
                             currRow = currRow + currRow2;
                             AUnknowns = [AUnknowns ; currRow]; %#ok<*AGROW>
                             bUnknowns = [bUnknowns ; 0];

                             % 4th derivative equality constraint
                             der4Traj = calcDerTraj(times, timeIndex, 4);
                             currRow = zeros(1, 8 * (numWP - 1));
                             currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der4Traj, 2)) ...
                             = der4Traj;
                             currRow2 = zeros(1, 8 * (numWP - 1));
                             der4Traj2 = -1 * der4Traj;
                             currRow2(1, 8*(i-2) + 1 : 8*(i-2) + size(der4Traj, 2)) ...
                             = der4Traj2;
                             currRow = currRow + currRow2;
                             AUnknowns = [AUnknowns ; currRow]; %#ok<*AGROW>
                             bUnknowns = [bUnknowns ; 0];     
                         %{
                         if acceleration is unknown, add two constraints: the acceleration
                         equality constraint and the 5th derivative constraint
                         %}             
                         elseif k ==2 && isnan(waypoints{1,i}(k, 1))
                             % acceleration equality constraint
                             der2Traj = calcDerTraj(times, timeIndex, 2);
                             currRow = zeros(1, 8 * (numWP - 1));
                             currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der2Traj, 2)) ...
                             = der2Traj;
                             currRow2 = zeros(1, 8 * (numWP - 1));
                             der2Traj2 = -1 * der2Traj;
                             currRow2(1, 8*(i-2) + 1 : 8*(i-2) + size(der2Traj, 2)) ...
                             = der2Traj2;
                             currRow = currRow + currRow2;
                             AUnknowns = [AUnknowns ; currRow]; %#ok<*AGROW>
                             bUnknowns = [bUnknowns ; 0];

                             % 5th derivative equality constraint
                             der5Traj = calcDerTraj(times, timeIndex, 5);
                             currRow = zeros(1, 8 * (numWP - 1));
                             currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der5Traj, 2)) ...
                             = der5Traj;
                             currRow2 = zeros(1, 8 * (numWP - 1));
                             der5Traj2 = -1 * der5Traj;
                             currRow2(1, 8*(i-2) + 1 : 8*(i-2) + size(der5Traj, 2)) ...
                             = der5Traj2;
                             currRow = currRow + currRow2;
                             AUnknowns = [AUnknowns ; currRow]; %#ok<*AGROW>
                             bUnknowns = [bUnknowns ; 0];
                         %{
                         if jerk is unknown, add two constraints: the jerk
                         equality constraint and the 6th derivative constraint
                         %}                 
                         elseif k == 3 && isnan(waypoints{1,i}(k, 1))
                             % jerk equality constraint
                             der3Traj = calcDerTraj(times, timeIndex, 3);
                             currRow = zeros(1, 8 * (numWP - 1));
                             currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der3Traj, 2)) ...
                             = der3Traj;
                             currRow2 = zeros(1, 8 * (numWP - 1));
                             der3Traj2 = -1 * der3Traj;
                             currRow2(1, 8*(i-2) + 1 : 8*(i-2) + size(der3Traj, 2)) ...
                             = der3Traj2;
                             currRow = currRow + currRow2;
                             AUnknowns = [AUnknowns ; currRow]; %#ok<*AGROW>
                             bUnknowns = [bUnknowns ; 0];

                             % 6th derivative equality constraint
                             der6Traj = calcDerTraj(times, timeIndex, 6);
                             currRow = zeros(1, 8 * (numWP - 1));
                             currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der6Traj, 2)) ...
                             = der6Traj;
                             currRow2 = zeros(1, 8 * (numWP - 1));
                             der6Traj2 = -1 * der6Traj;
                             currRow2(1, 8*(i-2) + 1 : 8*(i-2) + size(der6Traj, 2)) ...
                             = der6Traj2;
                             currRow = currRow + currRow2;
                             AUnknowns = [AUnknowns ; currRow]; %#ok<*AGROW>
                             bUnknowns = [bUnknowns ; 0];
                         end         
                    end
                end
            end
            b = [b ; bUnknowns];
            A = [A ; AUnknowns];
            x = A \ b;
            vals = {A, b, x};
        end
        
        function traj = getTrajectory(obj)
            valsCell = obj.generateTrajAndCoeff();
            traj = valsCell{1, size(vals, 2)};
        end
        
        function mat = getCoefficientMatrix(obj)
            valsCell = obj.generateTrajAndCoeff();
            mat = valsCell{1, 1};
        end
        
        function out = getOutputVec(obj)
            valsCell = obj.generateTrajAndCoeff();
            out = valsCell{1, size(vals, 2) - 1};
        end        
        
        % here, for all the plotting functions, x is trajectory polynomial
        function dummyVar = plotPosition(obj)
            x = obj.getTrajectory();
            times = obj.times;
            positions = obj.waypoints{1,1};
            numWP = size(positions, 1);
            numTraj = numWP - 1;
            sizeSol = size(x, 1);
            i = 1;
            timeIndex = 1;
            subplot(1,1,1);
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
            collectNames = cell(1, numTraj);
            for i = 1:numTraj
                text = strcat("trajectory ", num2str(i));
                collectNames{1,i} = text;
            end
            legend(collectNames);
            hold off
            dummyVar = [];
        end
        
        function dummyVar = plotVelocity(obj)
            x = obj.getTrajectory();
            times = obj.times;
            positions = obj.waypoints{1,1};
            numWP = size(positions, 1);
            numTraj = numWP - 1;
            sizeSol = size(x, 1);
            i = 1;
            timeIndex = 1;
            subplot(1,1,1);
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
            collectNames = cell(1, numTraj);
            for i = 1:numTraj
                text = strcat("trajectory ", num2str(i));
                collectNames{1,i} = text;
            end
            legend(collectNames);
            hold off
            dummyVar = [];
        end
        
        function dummyVar = plotAcceleration(obj)
            x = obj.getTrajectory();
            times = obj.times;
            positions = obj.waypoints{1,1};
            numWP = size(positions, 1);
            numTraj = numWP - 1;
            sizeSol = size(x, 1);
            i = 1;
            timeIndex = 1;
            subplot(1,1,1);
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
            collectNames = cell(1, numTraj);
            for i = 1:numTraj
                text = strcat("trajectory ", num2str(i));
                collectNames{1,i} = text;
            end
            legend(collectNames);
            hold off
            dummyVar = [];
        end
        
        function dummyVar = plotJerk(obj)
            x = obj.getTrajectory();
            times = obj.times;
            positions = obj.waypoints{1,1};
            numWP = size(positions, 1);
            numTraj = numWP - 1;
            sizeSol = size(x, 1);
            i = 1;
            timeIndex = 1;
            subplot(1,1,1);
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
            collectNames = cell(1, numTraj);
            for i = 1:numTraj
                text = strcat("trajectory ", num2str(i));
                collectNames{1,i} = text;
            end
            legend(collectNames);
            hold off
            dummyVar = [];
        end
        
        function dummyVar = plotAll(obj)
            x = obj.getTrajectory();
            times = obj.times;
            positions = obj.waypoints{1,1};
            numWP = size(positions, 1);
            numTraj = numWP - 1;
            sizeSol = size(x, 1);
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
            collectNames = cell(1, numTraj);
            for i = 1:numTraj
                text = strcat("trajectory ", num2str(i));
                collectNames{1,i} = text;
            end
            legend(collectNames);
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
            collectNames = cell(1, numTraj);
            for i = 1:numTraj
                text = strcat("trajectory ", num2str(i));
                collectNames{1,i} = text;
            end
            legend(collectNames);
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
            collectNames = cell(1, numTraj);
            for i = 1:numTraj
                text = strcat("trajectory ", num2str(i));
                collectNames{1,i} = text;
            end
            legend(collectNames);
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
            collectNames = cell(1, numTraj);
            for i = 1:numTraj
                text = strcat("trajectory ", num2str(i));
                collectNames{1,i} = text;
            end
            legend(collectNames);
            hold off
            
            dummyVar = [];
        end
        
        function derTraj = calcDerTraj(times, timeIndex, nthDer)
            % calculates derivatives as a row for trajectory
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
        
    end
    
end