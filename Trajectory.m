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
    
    properties (GetAccess = private)
        dimensions
        trajectory
    end
    
    methods (Abstract)
        % obj is like the 'this' keyword in Python or Java 
        dof = getDimensions(obj)
        traj = getTrajectory(obj)
        obj = setDimensions(obj)
        obj = setTrajectory(obj)
    end
end

