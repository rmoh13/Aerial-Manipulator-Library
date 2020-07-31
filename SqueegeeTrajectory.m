classdef SqueegeeTrajectory < Trajectory
    %SQUEEGEETRAJECTORY
    %{
    Generate trajectory to achieve orientation of squeegee, which is
    our aerial manipulator for the task of cleaning solar panels, on
    the quadrotor.
    
    %}
    
    properties (Access = private)
        Property1
    end
    
    methods
        function obj = SqueegeeTrajectory(inputArg1,inputArg2)
            %SQUEEGEETRAJECTORY Constructor
            %   Assume our Squeegee can move in all 3 Euler Angles 
            %   input is a matrix 
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function collectTrajs = generateTrajAndCoeff(obj, dim)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            collectTrajs = [];
        end
    end
end

