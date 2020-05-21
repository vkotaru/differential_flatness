classdef Flat2State < handle
% Class 
properties
    
    
end

%%
methods
    function obj = Flat2State(varargin)
        
        
    end   
end

methods (Static)
    % externally defined
    [ref] = computeQuadrotorMoment(aQ, daQ, d2aQ, varargin)
    [ref] = quadrotor(traj, mQ, J)
    [ref] = quadrotorLoad()
end
    
end