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
    
    properties (GetAccess = private)
        waypoints
        times
    end
    
    methods
        function obj = MultiSegmentTrajPlanner(waypoints,times)
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
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
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
    end
end

