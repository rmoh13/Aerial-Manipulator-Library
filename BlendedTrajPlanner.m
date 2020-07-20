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
            corresponding to them
            
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
            %}
            % do the first trajectory
            currentTrajObj = obj.cellArrayTrajs{1,1};
            dimensionsCurrentTraj = currentTrajObj.dimensions;
            currentdT = obj.dTs(1,1);
            endT = currentTrajObj.times(end, 1);
            t = [endT ; endT + currentdT];
            trajs = [];
            wp = {};
            for i = 1:dimensionsCurrentTraj
                trajs = [trajs ; currentTrajObj.getTrajectory(i)];
            end
            for i = 1:dimensionsCurrentTraj
                for j = 1:size(wp, 2)
                    
                end
            end
            
            
            % go in between
            for i = 2:size(obj.cellArrayTrajs)-1
                currentTrajObj = obj.cellArrayTrajs{1,i};
                dimensionsCurrentTraj = currentTrajObj.dimensions;
                
            end
            
            % do the last trajectory
            if size(obj.cellArrayTrajs) >= 2
                currentTrajObj = obj.cellArrayTrajs{1,end};
                dimensionsCurrentTraj = currentTrajObj.dimensions;
                currentdT = obj.dTs(end,1);
                startT = currentTrajObj.times(1, 1); % this is start of last trajectory/end of this one
                
                t = [endT ; endT + currentdT];
            end
            
        end
    end
end

