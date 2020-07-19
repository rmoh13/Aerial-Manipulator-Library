classdef (Abstract) Trajectory
    %TRAJECTORY - Abstract class as basis for different kinds of trajectory
    %objects.
    
    %{
    Properties - 
    # of degrees of freedom/dimensions 
    representation for the trajectory 
    (and this depends on whatever subclasses implement the interface 
    and how they represent the trajectory, for example, the multisegment planner 
    trajectory uses a polynomial, but we want to account for trajectories that 
    show the different angles like for the squeeguee trajectory and its rotation)
    
    Methods - getter, setter methods
    %}  
    
    properties (Access = protected)
        dimensions
        duration
    end
    
    methods
        function dof = get.dimensions(obj)
            dof = obj.dimensions;
        end
         
         function duration = get.duration(obj)
             duration = obj.duration;
         end
    end
    
    methods (Abstract)
        % obj is like the 'this' keyword in Python or Java 
        % no need to initialize getters here, they're implemented in subclasses
        traj = generateTrajAndCoeff(obj)
    end
end

