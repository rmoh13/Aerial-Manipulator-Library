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
            
            collectTrajs = {};
            trajs = [];
            %collectTrajs = {};
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
                 collectTrajs{end + 1} = genBlenTraj;
%                  for k = 1:dimensionsBeforeTraj
%                      trajs = [beforeTrajObj.getSpecificPositionTrajectory(endT); 
%                              genBlenTraj.getSpecificPositionTrajectory(t(1,1));
%                              afterTrajObj.getSpecificPositionTrajectory(afterTrajObj.times(1, 1))];
%                  end
                
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
            collectionTrajs = obj.generateTrajAndCoeff();
            generatedTrajs = collectionTrajs{1,1};
            for i = 1:size(collectionTrajs, 2)
                currentTrajectory = 
            end
            dummyVar = [];
        end
    end
end

