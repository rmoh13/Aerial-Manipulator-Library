classdef BlendedTrajPlanner < Trajectory
    %BLENDEDTRAJPLANNER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        cellArrayTrajs
        dTs
    end
    
    methods
        function obj = BlendedTrajPlanner(trajs, dTs)
            %BLENDEDTRAJPLANNER Constructor
            %{
            have a cell array in the constructor as given by 
            the user in the following order: a trajectory object 
            (it should be part of a class that extends the Trajectory 
            abstract class and thus has a getter method for a trajectory we can 
            use to eval positions, velocities, accelerations, jerks, etc.)
            and a column vector of dTs. And it must be alternating just in case 
            the user has multiple trajectory objects and sets of dTs for 
            corresponding to them. Assume there are at least 2 trajectory
            objects given by the user
            
            dT is the time difference that our generated trajectory
            operates between, so if dT = 5 then our generated trajectory
            goes from 
            %}
            obj.cellArrayTrajs = trajs;
            obj.dTs = dTs;
            obj.duration = obj.cellArrayTrajs{1, 1}.times(1,1) - ...
                           obj.cellArrayTrajs{1, size(obj.cellArrayTrajs, 2)}.times(end,1);
        end
        
        function cellArrayTrajs = get.cellArrayTrajs(obj)
            cellArrayTrajs = obj.cellArrayTrajs;
        end
        
        function dTs = get.dTs(obj)
            dTs = obj.dTs;
        end
        
        function collectTrajs = generateTrajAndCoeff(obj)
            % dim here is the total # of dimensions of the system
            %{
            loop through obj.cellArrayTrajs, evaluate the current traj for 
            the position, velocity, acceleration, jerk values at the first
            time point and evaluate it for the second time point which is
            going to be the first time point + the current dT and those are
            our ICs at those two points. We have a start and end point.
            
            output is collectTrajs where each x rows, where x is the 
            dimension, is a trajectory, and then we use this matrix to later plot
            the trajectories one by one (Taking into account their
            dimensions)
            %}
            
            tempCollectTrajs = {};
            %trajs = [];
            collectTrajs = {};
            wp = {};
            for i = 2:size(obj.cellArrayTrajs, 2)
                beforeTrajObj = obj.cellArrayTrajs{1, i-1};
                dimensionsBeforeTraj = beforeTrajObj.dimensions;
                currentdT = obj.dTs(i-1,1);
                endT = beforeTrajObj.times(end, 1);
                t = [0 ; currentdT];
                % we assume that all trajectories start from t = 0
                afterTrajObj = obj.cellArrayTrajs{1, i};
                % must be true: dimensionsAfterTraj == dimensionsBeforeTraj
                dimensionsAfterTraj = afterTrajObj.dimensions; %#ok<*NASGU>
                for j = 1:dimensionsBeforeTraj
                    j1 = beforeTrajObj.getPosition(endT);
                    j2 = afterTrajObj.getPosition(afterTrajObj.times(1, 1));
                    wp{j, 1} = [j1(j,1); 
                                j2(j,1)]; %#ok<*AGROW>
                    j1 = beforeTrajObj.getVelocity(endT);
                    j2 = beforeTrajObj.getAcceleration(endT);
                    j3 = beforeTrajObj.getJerk(endT);
                    wp{j, 2} = [j1(j,1); 
                                j2(j,1);
                                j3(j,1)];
                    j1 = afterTrajObj.getVelocity(afterTrajObj.times(1, 1));
                    j2 = afterTrajObj.getAcceleration(afterTrajObj.times(1, 1));
                    j3 = afterTrajObj.getJerk(afterTrajObj.times(1, 1));
                    wp{j, 3} = [j1(j,1); 
                                j2(j,1);
                                j3(j,1)];
                end
                % afterTrajObj is the one we're currenty on technically
                 genBlenTraj = MultiSegmentTrajPlanner(wp, t, dimensionsBeforeTraj);
                 tempCollectTrajs{end + 1} = beforeTrajObj;
                 tempCollectTrajs{end + 1} = genBlenTraj;
                 tempCollectTrajs{end + 1} = afterTrajObj;
%                  for k = 1:dimensionsBeforeTraj
%                      trajs = [beforeTrajObj.getSpecificPositionTrajectory(endT); 
%                              genBlenTraj.getSpecificPositionTrajectory(t(1,1));
%                              afterTrajObj.getSpecificPositionTrajectory(afterTrajObj.times(1, 1))];
%                  end
                
            end
            
            if size(tempCollectTrajs, 2) == 3
                collectTrajs{1,1} = tempCollectTrajs{1,1};
                collectTrajs{1,2} = tempCollectTrajs{1,2};
                collectTrajs{1,3} = tempCollectTrajs{1,3};
            elseif size(tempCollectTrajs, 2) > 3
                collectTrajs{1,1} = tempCollectTrajs{1,1};
                collectTrajs{1,2} = tempCollectTrajs{1,2};
                collectTrajs{1,3} = tempCollectTrajs{1,3};
                for i = 4:size(tempCollectTrajs, 2)
                    collectTrajs{end + 1} = tempCollectTrajs{1,i + 1};
                    collectTrajs{end + 1} = tempCollectTrajs{1,i + 2};
                    i = i + 3; %#ok<FXSET>
                end
            end
            
            
        end
        
        function dummyVar = plotBlendedPositionTraj(obj, dim)
            %{
            loop through the output of generateTrajAndCoeff(obj, dim),
            which is a cell array of the trajectories in order. Iterate
            over every x trajs where x is the dimension of any of these
            trajectory objects, plot.
 
            plot the first trajectory from its start time (should be 0)
            to its t_f (current global t_f). for the intermediate trajectories, 
            plot the values of the trajectory from its start time to its own t_f but
            against the values of t that go from current global t_f plus 
            the duration of the current and the
            %}
            %collectionTrajs = obj.generateTrajAndCoeff();
            %generatedTrajs = collectionTrajs{1,1};
            collectionTrajs = obj.generateTrajAndCoeff();
            if dim == 1
                subplot(1,1,1);
                hold on
                %for i = 1:size(collectionTrajs, 2) + size(obj.cellArrayTrajs, 2)
                for i = 1:size(collectionTrajs, 2)
                    currentTrajectory = collectionTrajs{1,i};
                    x = currentTrajectory.getTrajectory(dim);
                    times = currentTrajectory.times;
                    positions = currentTrajectory.waypoints{dim,1};
                    numWP = size(positions, 1);
                    numTraj = numWP - 1;
                    sizeSol = size(x, 1);
                    j = 1;
                    timeIndex = 1;
                    while j < sizeSol
                        posTraj = x(j:j+7).';
                        t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                        y = posTraj;
                        plot(t, polyval(y,t))
                        title("position vs time")
                        xlabel("t")
                        ylabel("x(t)")
                        j = j + 8;
                        timeIndex = timeIndex + 1;
                    end
                end
                collectNames = cell(1, numTraj);
                for k = 1:numTraj
                    text = strcat("trajectory ", num2str(k));
                    collectNames{1,k} = text;
                end
                legend(collectNames);
                hold off
                dummyVar = [];
            elseif dim == 2
                subplot(1,1,1);
                hold on
                view(3);
                for i = 1:size(collectionTrajs, 2)
                    currentTrajectory = collectionTrajs{1,i};
                    x = currentTrajectory.getTrajectory(1);
                    y = currentTrajectory.getTrajectory(2);
                    times = currentTrajectory.times;
                    xpositions = currentTrajectory.waypoints{1,1};
                    ypositions = currentTrajectory.waypoints{2,1};
                    numWPx = size(xpositions, 1);
                    % these two should be equal
                    numWPy = size(ypositions, 1);
                    numXTraj = numWPx - 1;
                    % these two should be equal
                    numYTraj = numWPy - 1; %#ok<*NASGU>
                    sizeXSol = size(x, 1);
                    % these two should be equal
                    sizeYSol = size(y, 1);
                    j = 1;
                    timeIndex = 1;
                    while j < sizeXSol
                        xposTraj = x(j:j+7).';
                        yposTraj = y(j:j+7).';
                        t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                        xtraj = xposTraj;
                        ytraj = yposTraj;
                        plot(polyval(xtraj,t), polyval(ytraj,t))
                        title("position vs time")
                        xlabel("x(t)")
                        ylabel("y(t)")
                        j = j + 8;
                        timeIndex = timeIndex + 1;
                    end
                end
                collectNames = cell(1, numXTraj);
                for k = 1:numXTraj
                    text = strcat("trajectory ", num2str(k));
                    collectNames{1,k} = text;
                end
                legend(collectNames);
                hold off
                dummyVar = [];
            elseif dim == 3
                subplot(1,1,1);
                hold on
                view(3);
                for i = 1:size(collectionTrajs, 2)
                    currentTrajectory = collectionTrajs{1,i};
                    x = currentTrajectory.getTrajectory(1);
                    y = currentTrajectory.getTrajectory(2);
                    z = currentTrajectory.getTrajectory(3);
                    times = currentTrajectory.times;
                    xpositions = currentTrajectory.waypoints{1,1};
                    ypositions = currentTrajectory.waypoints{2,1};
                    zpositions = currentTrajectory.waypoints{3,1};
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
                    j = 1;
                    timeIndex = 1;
                    while j < sizeXSol
                        xposTraj = x(j:j+7).';
                        yposTraj = y(j:j+7).';
                        zposTraj = z(j:j+7).';
                        t = linspace(times(timeIndex,1),times(timeIndex + 1, 1));
                        xtraj = xposTraj;
                        ytraj = yposTraj;
                        ztraj = zposTraj;
                        plot3(polyval(xtraj,t), polyval(ytraj,t), polyval(ztraj,t))
                        title("position vs time")
                        xlabel("x(t)")
                        ylabel("y(t)")
                        zlabel("z(t)")
                        j = j + 8;
                        timeIndex = timeIndex + 1;
                    end
                end
                collectNames = cell(1, numXTraj);
                for k = 1:numXTraj
                    text = strcat("trajectory ", num2str(k));
                    collectNames{1,k} = text;
                end
                legend(collectNames);
                hold off
                dummyVar = [];
            end
        end
    end
end

