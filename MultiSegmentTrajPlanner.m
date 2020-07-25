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
    
    properties
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
            the first row in waypoints is for the 1st dimension
            the second row in waypoints is for the 2nd dimension,
            etc.
            
            EXAMPLE 1D 3 WAYPOINTS CELL ARRAY:
            {[0;5;10],[0;0;0],[10;20;30],[0;0;0]}
            
            
            EXAMPLE 2D 3 WAYPOINTS CELL ARRAY:
            {[0;5;10],[0;0;0],[10;20;30],[0;0;0];
            [0;7;11],[0;0;0],[1;2;3],[0;0;0]}
            
            EXAMPLE 3D 3 WAYPOINTS CELL ARRAY:
            {[0;5;10],[0;0;0],[10;20;30],[0;0;0];
            [0;7;11],[0;0;0],[1;2;3],[0;0;0];
            [0;8;100], [0;0;0], [4;5;6], [0;0;0]}
            
            

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
            % assume duration is from start t to end t that user specifies
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
        
        function derTraj = calcDerTraj(obj, times, timeIndex, nthDer) %#ok<INUSL>
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
        
        function vals = generateTrajAndCoeff(obj, dim)
            % in Ax = b, vals is x which is the constraints column vector
            % make waypoints rows, assume each position is given
            waypoints = obj.waypoints; %#ok<*PROPLC,*PROP>
            times = obj.times;
            positions = waypoints{dim,1};
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
            bVals = [waypoints{dim,2}];
            AValsFirstRow = [];
            % do the first row of the corresponding value constraints in A
            for i = 1:size(waypoints{dim, 2}, 1)
                if i == 1 && ~isnan(waypoints{dim,2}(i, 1))
                    AValsFirstRow = [AValsFirstRow ; obj.calcDerTraj(times, 1, 1) , ...  
                    zeros(1, size(A, 2) - size(obj.calcDerTraj(times, 1, 1), 2))]; 
                elseif i == 2 && ~isnan(waypoints{dim,2}(i, 1))
                    AValsFirstRow = [AValsFirstRow ; obj.calcDerTraj(times, 1, 2) , ...
                    zeros(1, size(A, 2) - size(obj.calcDerTraj(times, 1, 2), 2))]; 
                elseif i == 3 && ~isnan(waypoints{dim,2}(i, 1))
                    AValsFirstRow = [AValsFirstRow ; obj.calcDerTraj(times, 1, 3) , ... 
                    zeros(1, size(A, 2) - size(obj.calcDerTraj(times, 1, 3), 2))]; 
                end
            end
            A = [A ; AValsFirstRow];
        
            % now go in between
            j = 3;
            while j < size(waypoints, 2)
                currentCol = [];
                for k = 1:size(waypoints{dim,j}, 1)
                    if ~isnan(waypoints{dim,j}(k, 1))
                        currentCol = [currentCol ; waypoints{dim,j}(k, 1)];
                    end
                end
                bVals = [bVals ; currentCol; currentCol]; 
                j = j + 1;
            end
            bVals = [bVals ; waypoints{dim, size(waypoints, 2)}];
            b = [b ; bVals];
            AInBetween = [];
            for i = 3:size(waypoints, 2) - 1
                timeIndex = i - 1;
                for k = 1:size(waypoints{dim, i}, 1)
                     if k == 1 && ~isnan(waypoints{dim,i}(k, 1))
                         der1Traj = obj.calcDerTraj(times, timeIndex, 1);
                         currRow = zeros(1, 8 * (numWP - 1));
                         currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der1Traj, 2)) ...
                         = der1Traj;
                         AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
                     elseif k == 2 && ~isnan(waypoints{dim,i}(k, 1))
                         der2Traj = obj.calcDerTraj(times, timeIndex, 2);
                         currRow = zeros(1, 8 * (numWP - 1));
                         currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der2Traj, 2)) ...
                         = der2Traj;
                         AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
                     elseif k == 3 && ~isnan(waypoints{dim,i}(k, 1))
                         der3Traj = obj.calcDerTraj(times, timeIndex, 3);
                         currRow = zeros(1, 8 * (numWP - 1));
                         currRow(1, 8*(i-3) + 1 : 8*(i-3) + size(der3Traj, 2)) ...
                         = der3Traj;
                         AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
                     end         
                end
                for k = 1:size(waypoints{dim, i}, 1)
                     if k == 1 && ~isnan(waypoints{dim,i}(k, 1))
                         der1Traj = obj.calcDerTraj(times, timeIndex, 1);
                         currRow = zeros(1, 8 * (numWP - 1));
                         currRow(1, 8*(i-2) + 1 : 8*(i-2) + size(der1Traj, 2)) ...
                         = der1Traj;
                         AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
                     elseif k == 2 && ~isnan(waypoints{dim,i}(k, 1))
                         der2Traj = obj.calcDerTraj(times, timeIndex, 2);
                         currRow = zeros(1, 8 * (numWP - 1));
                         currRow(1, 8*(i-2) + 1 : 8*(i-2) + size(der2Traj, 2)) ...
                         = der2Traj;
                         AInBetween = [AInBetween ; currRow]; %#ok<*AGROW>
                     elseif k == 3 && ~isnan(waypoints{dim,i}(k, 1))
                         der3Traj = obj.calcDerTraj(times, timeIndex, 3);
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
            for k = 1:size(waypoints{dim, size(waypoints, 2)}, 1)
                 if k == 1 && ~isnan(waypoints{dim,size(waypoints, 2)}(k, 1))
                     der1Traj = obj.calcDerTraj(times, lastWNaNndex, 1);
                     currRow = zeros(1, 8 * (numWP - 1));
                     currRow(1, 8*(lastWNaNndex-2) + 1 : 8*(lastWNaNndex-2) + size(der1Traj, 2)) ...
                     = der1Traj;
                     AValsLastRow = [AValsLastRow ; currRow]; %#ok<*AGROW>
                 elseif k == 2 && ~isnan(waypoints{dim,size(waypoints, 2)}(k, 1))
                     der2Traj = obj.calcDerTraj(times, lastWNaNndex, 2);
                     currRow = zeros(1, 8 * (numWP - 1));
                     currRow(1, 8*(lastWNaNndex-2) + 1 : 8*(lastWNaNndex-2) + size(der2Traj, 2)) ...
                     = der2Traj;
                     AValsLastRow = [AValsLastRow ; currRow]; %#ok<*AGROW>
                 elseif k == 3 && ~isnan(waypoints{dim,size(waypoints, 2)}(k, 1))
                     der3Traj = obj.calcDerTraj(times, lastWNaNndex, 3);
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
                indicesUnknowns = find(isnan(waypoints{dim, i}));
                if size(indicesUnknowns, 1) ~= 0
                    for k = 1:size(waypoints{dim, i}, 1)
                        %currentIndexUnknown = indicesUnknowns(1, k);
                        %waypoints{1,i}(currentIndexUnknown, 1);
                         %{
                         if velocity is unknown, add two constraints: the velocity
                         equality constraint and the 4th derivative constraint
                         %}
                         if k == 1 && isnan(waypoints{dim,i}(k, 1))
                             % velocity equality constraint
                             der1Traj = obj.calcDerTraj(times, timeIndex, 1);
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
                             der4Traj = obj.calcDerTraj(times, timeIndex, 4);
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
                         elseif k ==2 && isnan(waypoints{dim,i}(k, 1))
                             % acceleration equality constraint
                             der2Traj = obj.calcDerTraj(times, timeIndex, 2);
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
                             der5Traj = obj.calcDerTraj(times, timeIndex, 5);
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
                         elseif k == 3 && isnan(waypoints{dim,i}(k, 1))
                             % jerk equality constraint
                             der3Traj = obj.calcDerTraj(times, timeIndex, 3);
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
                             der6Traj = obj.calcDerTraj(times, timeIndex, 6);
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
        
        function traj = getTrajectory(obj, dim)
            valsCell = obj.generateTrajAndCoeff(dim);
            traj = valsCell{1, size(valsCell, 2)};
        end
        
        function mat = getCoefficientMatrix(obj, dim)
            valsCell = obj.generateTrajAndCoeff(dim);
            mat = valsCell{1, 1};
        end
        
        function out = getOutputVec(obj, dim)
            valsCell = obj.generateTrajAndCoeff(dim);
            out = valsCell{1, size(valsCell, 2) - 1};
        end
        
        function traj = getSpecificPositionTrajectory(obj, time)
            if obj.dimensions == 1
                x = obj.getTrajectory(1);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    posTraj = x(i:i+7).';
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = posTraj;
                    end
                end
            elseif obj.dimensions == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = [xposTraj ; yposTraj];
                    end
                end
            elseif obj.dimensions == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    zposTraj = z(i:i+7).';
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = [xposTraj ; yposTraj ; zposTraj];
                    end
                end
            end
        end
        
        function traj = getSpecificVelocityTrajectory(obj, time)
            if obj.dimensions == 1
                x = obj.getTrajectory(1);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    posTraj = x(i:i+7).';
                    velTraj = polyder(posTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = velTraj;
                    end
                end
            elseif obj.dimensions == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = [xvelTraj ; yvelTraj];
                    end
                end
            elseif obj.dimensions == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    zposTraj = z(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    zvelTraj = polyder(zposTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = [xvelTraj ; yvelTraj ; zvelTraj];
                    end
                end
            end
        end
        
        function traj = getSpecificAccelerationTrajectory(obj, time)
            if obj.dimensions == 1
                x = obj.getTrajectory(1);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    posTraj = x(i:i+7).';
                    velTraj = polyder(posTraj);
                    accelTraj = polyder(velTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = accelTraj;
                    end
                end
            elseif obj.dimensions == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = [xaccelTraj ; yaccelTraj];
                    end
                end
            elseif obj.dimensions == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    zposTraj = z(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    zvelTraj = polyder(zposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    zaccelTraj = polyder(zvelTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = [xaccelTraj ; yaccelTraj ; zaccelTraj];
                    end
                end
            end
        end
        
        function traj = getSpecificJerkTrajectory(obj, time)
            if obj.dimensions == 1
                x = obj.getTrajectory(1);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    posTraj = x(i:i+7).';
                    velTraj = polyder(posTraj);
                    accelTraj = polyder(velTraj);
                    jerkTraj = polyder(accelTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = jerkTraj;
                    end
                end
            elseif obj.dimensions == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    xjerkTraj = polyder(xaccelTraj);
                    yjerkTraj = polyder(yaccelTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = [xjerkTraj ; yjerkTraj];
                    end
                end
            elseif obj.dimensions == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    zposTraj = z(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    zvelTraj = polyder(zposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    zaccelTraj = polyder(zvelTraj);
                    xjerkTraj = polyder(xaccelTraj);
                    yjerkTraj = polyder(yaccelTraj);
                    zjerkTraj = polyder(zaccelTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        traj = [xjerkTraj ; yjerkTraj ; zjerkTraj];
                    end
                end
            end
        end
        
        function vec = getPosition(obj, time)
            if obj.dimensions == 1
                x = obj.getTrajectory(1);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    posTraj = x(i:i+7).';
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = polyval(posTraj, time);
                    end
                end
            elseif obj.dimensions == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = [polyval(xposTraj, time) ; polyval(yposTraj, time)];
                    end
                end
            elseif obj.dimensions == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    zposTraj = z(i:i+7).';
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = [polyval(xposTraj, time) ; polyval(yposTraj, time) ; polyval(zposTraj, time)];
                    end
                end
            end
        end
        
        function vec = getVelocity(obj, time)
            if obj.dimensions == 1
                x = obj.getTrajectory(1);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    posTraj = x(i:i+7).';
                    velTraj = polyder(posTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = polyval(velTraj, time);
                    end
                end
            elseif obj.dimensions == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = [polyval(xvelTraj, time) ; polyval(yvelTraj, time)];
                    end
                end
            elseif obj.dimensions == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    zposTraj = z(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    zvelTraj = polyder(zposTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = [polyval(xvelTraj, time) ; polyval(yvelTraj, time) ; polyval(zvelTraj, time)];
                    end
                end
            end
        end
        
        function vec = getAcceleration(obj, time)
            if obj.dimensions == 1
                x = obj.getTrajectory(1);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    posTraj = x(i:i+7).';
                    velTraj = polyder(posTraj);
                    accelTraj = polyder(velTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = polyval(accelTraj, time);
                    end
                end
            elseif obj.dimensions == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = [polyval(xaccelTraj, time) ; polyval(yaccelTraj, time)];
                    end
                end
            elseif obj.dimensions == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    zposTraj = z(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    zvelTraj = polyder(zposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    zaccelTraj = polyder(zvelTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = [polyval(xaccelTraj, time) ; polyval(yaccelTraj, time) ; polyval(zaccelTraj, time)];
                    end
                end
            end
        end
        
        function vec = getJerk(obj, time)
            if obj.dimensions == 1
                x = obj.getTrajectory(1);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    posTraj = x(i:i+7).';
                    velTraj = polyder(posTraj);
                    accelTraj = polyder(velTraj);
                    jerkTraj = polyder(accelTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = polyval(jerkTraj, time);
                    end
                end
            elseif obj.dimensions == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    xjerkTraj = polyder(xaccelTraj);
                    yjerkTraj = polyder(yaccelTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = [polyval(xjerkTraj, time) ; polyval(yjerkTraj, time)];
                    end
                end
            elseif obj.dimensions == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                sizeSol = size(x, 1);
                timeIndex = 1;
                i = 1;
                while i < sizeSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    zposTraj = z(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yvelTraj = polyder(yposTraj);
                    zvelTraj = polyder(zposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    zaccelTraj = polyder(zvelTraj);
                    xjerkTraj = polyder(xaccelTraj);
                    yjerkTraj = polyder(yaccelTraj);
                    zjerkTraj = polyder(zaccelTraj);
                    timeInterval = linspace(obj.times(timeIndex,1),obj.times(timeIndex + 1, 1));
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                    if time >= timeInterval(1,1) && time <= timeInterval(1, end)
                        vec = [polyval(xjerkTraj, time) ; polyval(yjerkTraj, time) ; polyval(zjerkTraj, time)];
                    end
                end
            end
        end
        
        % here, for all the plotting functions, x is trajectory polynomial
        function dummyVar = plotPosition(obj, dim)
            if dim == 1
                x = obj.getTrajectory(dim);
                times = obj.times;
                positions = obj.waypoints{dim,1};
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
            elseif dim == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                times = obj.times;
                xpositions = obj.waypoints{1,1};
                ypositions = obj.waypoints{2,1};
                numWPx = size(xpositions, 1);
                % these two should be equal
                numWPy = size(ypositions, 1);
                numXTraj = numWPx - 1;
                % these two should be equal
                numYTraj = numWPy - 1; %#ok<*NASGU>
                sizeXSol = size(x, 1);
                % these two should be equal
                sizeYSol = size(y, 1);
                i = 1;
                timeIndex = 1;
                subplot(1,1,1);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    xtraj = xposTraj;
                    ytraj = yposTraj;
                    plot(polyval(xtraj,t), polyval(ytraj,t))
                    title("position vs time")
                    xlabel("x(t)")
                    ylabel("y(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off
                dummyVar = [];
            elseif dim == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                times = obj.times;
                xpositions = obj.waypoints{1,1};
                ypositions = obj.waypoints{2,1};
                zpositions = obj.waypoints{3,1};
                numWPx = size(xpositions, 1);
                % these three should be equal
                numWPy = size(ypositions, 1);
                numWPz = size(zpositions, 1);
                numXTraj = numWPx - 1;
                % these three should be equal
                numYTraj = numWPy - 1; %#ok<*NASGU>
                numZTraj = numWPz - 1; %#ok<*NASGU>
                sizeXSol = size(x, 1);
                % these three should be equal
                sizeYSol = size(y, 1);
                sizeZSol = size(z, 1);
                i = 1;
                timeIndex = 1;
                subplot(1,1,1);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    zposTraj = z(i:i+7).';
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    xtraj = xposTraj;
                    ytraj = yposTraj;
                    ztraj = zposTraj;
                    plot3(polyval(xtraj,t), polyval(ytraj,t), polyval(ztraj,t))
                    title("position vs time")
                    xlabel("x(t)")
                    ylabel("y(t)")
                    zlabel("z(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off
                dummyVar = [];
            end
        end
        
        function dummyVar = plotVelocity(obj, dim)
            if dim == 1
                x = obj.getTrajectory(dim);
                times = obj.times;
                positions = obj.waypoints{dim,1};
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
            elseif dim == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                times = obj.times;
                xpositions = obj.waypoints{1,1};
                ypositions = obj.waypoints{2,1};
                numWPx = size(xpositions, 1);
                % these two should be equal
                numWPy = size(ypositions, 1);
                numXTraj = numWPx - 1;
                % these two should be equal
                numYTraj = numWPy - 1; %#ok<*NASGU>
                sizeXSol = size(x, 1);
                % these two should be equal
                sizeYSol = size(y, 1);
                i = 1;
                timeIndex = 1;
                subplot(1,1,1);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot(polyval(xvelTraj,t), polyval(yvelTraj,t))
                    title("velocity vs time")
                    xlabel("x'(t)")
                    ylabel("y'(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off
                dummyVar = [];
            elseif dim == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                times = obj.times;
                xpositions = obj.waypoints{1,1};
                ypositions = obj.waypoints{2,1};
                zpositions = obj.waypoints{3,1};
                numWPx = size(xpositions, 1);
                % these three should be equal
                numWPy = size(ypositions, 1);
                numWPz = size(zpositions, 1);
                numXTraj = numWPx - 1;
                % these three should be equal
                numYTraj = numWPy - 1; %#ok<*NASGU>
                numZTraj = numWPz - 1; %#ok<*NASGU>
                sizeXSol = size(x, 1);
                % these three should be equal
                sizeYSol = size(y, 1);
                sizeZSol = size(z, 1);
                i = 1;
                timeIndex = 1;
                subplot(1,1,1);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    zposTraj = z(i:i+7).';
                    zvelTraj = polyder(zposTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot3(polyval(xvelTraj,t), polyval(yvelTraj,t), polyval(zvelTraj,t))
                    title("velocity vs time")
                    xlabel("x'(t)")
                    ylabel("y'(t)")
                    zlabel("z'(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off
                dummyVar = [];
            end
        end
        
        function dummyVar = plotAcceleration(obj, dim)
            if dim == 1
                x = obj.getTrajectory(dim);
                times = obj.times;
                positions = obj.waypoints{dim,1};
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
            elseif dim == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                times = obj.times;
                xpositions = obj.waypoints{1,1};
                ypositions = obj.waypoints{2,1};
                numWPx = size(xpositions, 1);
                % these two should be equal
                numWPy = size(ypositions, 1);
                numXTraj = numWPx - 1;
                % these two should be equal
                numYTraj = numWPy - 1; %#ok<*NASGU>
                sizeXSol = size(x, 1);
                % these two should be equal
                sizeYSol = size(y, 1);
                i = 1;
                timeIndex = 1;
                subplot(1,1,1);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot(polyval(xaccelTraj,t), polyval(yaccelTraj,t))
                    title("acceleration vs time")
                    xlabel("x''(t)")
                    ylabel("y''(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off
                dummyVar = [];
            elseif dim == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                times = obj.times;
                xpositions = obj.waypoints{1,1};
                ypositions = obj.waypoints{2,1};
                zpositions = obj.waypoints{3,1};
                numWPx = size(xpositions, 1);
                % these three should be equal
                numWPy = size(ypositions, 1);
                numWPz = size(zpositions, 1);
                numXTraj = numWPx - 1;
                % these three should be equal
                numYTraj = numWPy - 1; %#ok<*NASGU>
                numZTraj = numWPz - 1; %#ok<*NASGU>
                sizeXSol = size(x, 1);
                % these three should be equal
                sizeYSol = size(y, 1);
                sizeZSol = size(z, 1);
                i = 1;
                timeIndex = 1;
                subplot(1,1,1);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    zposTraj = z(i:i+7).';
                    zvelTraj = polyder(zposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    zaccelTraj = polyder(zvelTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot3(polyval(xaccelTraj,t), polyval(yaccelTraj,t), polyval(zaccelTraj,t))
                    title("acceleration vs time")
                    xlabel("x''(t)")
                    ylabel("y''(t)")
                    zlabel("z''(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off
                dummyVar = [];
            end
        end
        
        function dummyVar = plotJerk(obj, dim)
            if dim == 1
                x = obj.getTrajectory(dim);
                times = obj.times;
                positions = obj.waypoints{dim,1};
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
            elseif dim == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                times = obj.times;
                xpositions = obj.waypoints{1,1};
                ypositions = obj.waypoints{2,1};
                numWPx = size(xpositions, 1);
                % these two should be equal
                numWPy = size(ypositions, 1);
                numXTraj = numWPx - 1;
                % these two should be equal
                numYTraj = numWPy - 1; %#ok<*NASGU>
                sizeXSol = size(x, 1);
                % these two should be equal
                sizeYSol = size(y, 1);
                i = 1;
                timeIndex = 1;
                subplot(1,1,1);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    xjerkTraj = polyder(xaccelTraj);
                    yjerkTraj = polyder(yaccelTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot(polyval(xjerkTraj,t), polyval(yjerkTraj,t))
                    title("jerk vs time")
                    xlabel("x'''(t)")
                    ylabel("y'''(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off
                dummyVar = [];
            elseif dim == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                times = obj.times;
                xpositions = obj.waypoints{1,1};
                ypositions = obj.waypoints{2,1};
                zpositions = obj.waypoints{3,1};
                numWPx = size(xpositions, 1);
                % these three should be equal
                numWPy = size(ypositions, 1);
                numWPz = size(zpositions, 1);
                numXTraj = numWPx - 1;
                % these three should be equal
                numYTraj = numWPy - 1; %#ok<*NASGU>
                numZTraj = numWPz - 1; %#ok<*NASGU>
                sizeXSol = size(x, 1);
                % these three should be equal
                sizeYSol = size(y, 1);
                sizeZSol = size(z, 1);
                i = 1;
                timeIndex = 1;
                subplot(1,1,1);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    zposTraj = z(i:i+7).';
                    zvelTraj = polyder(zposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    zaccelTraj = polyder(zvelTraj);
                    xjerkTraj = polyder(xaccelTraj);
                    yjerkTraj = polyder(yaccelTraj);
                    zjerkTraj = polyder(zaccelTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot3(polyval(xjerkTraj,t), polyval(yjerkTraj,t), polyval(zjerkTraj,t))
                    title("jerk vs time")
                    xlabel("x'''(t)")
                    ylabel("y'''(t)")
                    zlabel("z'''(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off
                dummyVar = [];
            end
        end
        
        % this is basically for 1D only. dim isn't even necessary here
        function dummyVar = plot1DimAll(obj, dim)
            x = obj.getTrajectory(dim);
            times = obj.times;
            positions = obj.waypoints{dim,1};
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
        
        function dummyVar = plotMultiDimAll(obj, dim)
            if dim == 2
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                times = obj.times;
                xpositions = obj.waypoints{1,1};
                ypositions = obj.waypoints{2,1};
                numWPx = size(xpositions, 1);
                % these two should be equal
                numWPy = size(ypositions, 1);
                numXTraj = numWPx - 1;
                % these two should be equal
                numYTraj = numWPy - 1; %#ok<*NASGU>
                sizeXSol = size(x, 1);
                % these two should be equal
                sizeYSol = size(y, 1);
                i = 1;
                timeIndex = 1;
                subplot(4,1,1);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    xtraj = xposTraj;
                    ytraj = yposTraj;
                    plot(polyval(xtraj,t), polyval(ytraj,t))
                    title("position vs time")
                    xlabel("x(t)")
                    ylabel("y(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off

                i = 1;
                timeIndex = 1;
                subplot(4,1,2);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot(polyval(xvelTraj,t), polyval(yvelTraj,t))
                    title("velocity vs time")
                    xlabel("x'(t)")
                    ylabel("y'(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off

                i = 1;
                timeIndex = 1;
                subplot(4,1,3);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot(polyval(xaccelTraj,t), polyval(yaccelTraj,t))
                    title("acceleration vs time")
                    xlabel("x''(t)")
                    ylabel("y''(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off

                i = 1;
                timeIndex = 1;
                subplot(4,1,4);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    xjerkTraj = polyder(xaccelTraj);
                    yjerkTraj = polyder(yaccelTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot(polyval(xjerkTraj,t), polyval(yjerkTraj,t))
                    title("jerk vs time")
                    xlabel("x'''(t)")
                    ylabel("y'''(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off
            elseif dim == 3
                x = obj.getTrajectory(1);
                y = obj.getTrajectory(2);
                z = obj.getTrajectory(3);
                times = obj.times;
                xpositions = obj.waypoints{1,1};
                ypositions = obj.waypoints{2,1};
                zpositions = obj.waypoints{3,1};
                numWPx = size(xpositions, 1);
                % these three should be equal
                numWPy = size(ypositions, 1);
                numWPz = size(zpositions, 1);
                numXTraj = numWPx - 1;
                % these three should be equal
                numYTraj = numWPy - 1; %#ok<*NASGU>
                numZTraj = numWPz - 1; %#ok<*NASGU>
                sizeXSol = size(x, 1);
                % these three should be equal
                sizeYSol = size(y, 1);
                sizeZSol = size(z, 1);
                i = 1;
                timeIndex = 1;
                subplot(4,1,1);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    yposTraj = y(i:i+7).';
                    zposTraj = z(i:i+7).';
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    xtraj = xposTraj;
                    ytraj = yposTraj;
                    ztraj = zposTraj;
                    plot3(polyval(xtraj,t), polyval(ytraj,t), polyval(ztraj,t))
                    title("position vs time")
                    xlabel("x(t)")
                    ylabel("y(t)")
                    zlabel("z(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off

                i = 1;
                timeIndex = 1;
                subplot(4,1,2);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    zposTraj = z(i:i+7).';
                    zvelTraj = polyder(zposTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot3(polyval(xvelTraj,t), polyval(yvelTraj,t), polyval(zvelTraj,t))
                    title("velocity vs time")
                    xlabel("x'(t)")
                    ylabel("y'(t)")
                    zlabel("z'(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off

                i = 1;
                timeIndex = 1;
                subplot(4,1,3);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    zposTraj = z(i:i+7).';
                    zvelTraj = polyder(zposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    zaccelTraj = polyder(zvelTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot3(polyval(xaccelTraj,t), polyval(yaccelTraj,t), polyval(zaccelTraj,t))
                    title("acceleration vs time")
                    xlabel("x''(t)")
                    ylabel("y''(t)")
                    zlabel("z''(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off

                i = 1;
                timeIndex = 1;
                subplot(4,1,4);
                hold on
                view(3);
                while i < sizeXSol
                    xposTraj = x(i:i+7).';
                    xvelTraj = polyder(xposTraj);
                    yposTraj = y(i:i+7).';
                    yvelTraj = polyder(yposTraj);
                    zposTraj = z(i:i+7).';
                    zvelTraj = polyder(zposTraj);
                    xaccelTraj = polyder(xvelTraj);
                    yaccelTraj = polyder(yvelTraj);
                    zaccelTraj = polyder(zvelTraj);
                    xjerkTraj = polyder(xaccelTraj);
                    yjerkTraj = polyder(yaccelTraj);
                    zjerkTraj = polyder(zaccelTraj);
                    t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                    plot3(polyval(xjerkTraj,t), polyval(yjerkTraj,t), polyval(zjerkTraj,t))
                    title("jerk vs time")
                    xlabel("x'''(t)")
                    ylabel("y'''(t)")
                    zlabel("z'''(t)")
                    i = i + 8;
                    timeIndex = timeIndex + 1;
                end
                collectNames = cell(1, numXTraj);
                for i = 1:numXTraj
                    text = strcat("trajectory ", num2str(i));
                    collectNames{1,i} = text;
                end
                legend(collectNames);
                hold off
            end
            dummyVar = [];
        end
        
    end
    
end