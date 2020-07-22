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
        end
        
        function cellArrayTrajs = get.cellArrayTrajs(obj)
            cellArrayTrajs = obj.cellArrayTrajs;
        end
        
        function dTs = get.dTs(obj)
            dTs = obj.dTs;
        end
        
        function collectTrajs = generateTrajAndCoeff(obj, dim)
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
            
            for i = 2:size(obj.cellArrayTrajs, 2)-1
                collectTrajs = [];
                beforeTrajObj = obj.cellArrayTrajs{i-1,1};
                dimensionsCurrentTraj = currentTrajObj.dimensions;
                currentdT = obj.dTs(1,1);
                endT = currentTrajObj.times(end, 1);
                t = [endT ; endT + currentdT];
                trajs = [];
                wp = {};
                for i = 1:dimensionsCurrentTraj
                    trajectoryRow = currentTrajObj.getTrajectory(i);
                    trajectoryRow = trajectoryRow(end - 7:end, 1).';
                    trajs = [trajs ; trajectoryRow] %#ok<*AGROW>
                end
                afterTrajs = [];
                afterTrajObj = obj.cellArrayTrajs{1,2};
                for i = 1:dimensionsCurrentTraj
                    afterTrajs = [afterTrajs ; afterTrajObj.getTrajectory(i)]; %#ok<*AGROW>
                end
                for i = 1:dimensionsCurrentTraj
                    currentTraj = trajs(i, 1:end);
                    collectTrajs = [collectTrajs ; currentTraj];
                    if size(obj.cellArrayTrajs) > 1
                        currentNextTraj = nextTrajs(i, 1:end);
                    end

                    wp{i, 1} = [polyval(currentTraj, t(1,1)) ; polyval(currentNextTraj, t(2,1))];

                    wp{i, 2} = ...
                        [polyval(polyder(currentTraj), t(1,1)) ; 
                        polyval(polyder(polyder(currentTraj)), t(1,1)) ;
                        polyval(polyder(polyder(polyder(currentTraj))), t(1,1))];                

                    size(obj.cellArrayTrajs) > 1
                    wp{i, 3} = ...
                        [polyval(polyder(currentNextTraj), t(2,1)) ; 
                        polyval(polyder(polyder(currentNextTraj)), t(2,1)) ;
                        polyval(polyder(polyder(polyder(currentNextTraj))), t(2,1))];
                end
                celldisp(wp)
                t
                genBlenTraj = MultiSegmentTrajPlanner(wp, t, dimensionsCurrentTraj);
                genBlenTrajs = [];
                for i = 1:dimensionsCurrentTraj
                    genBlenTrajs = [genBlenTrajs ; genBlenTraj.getTrajectory(i)];
                    collectTrajs = [collectTrajs ; genBlenTrajs];
                end
            
            end
            
        end
    end
end

