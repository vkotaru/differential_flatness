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
    [ref] = quadrotor(traj, mQ, J)
    [ref] = computeQuadrotorMoment(aQ, daQ, d2aQ, varargin)
    [ref] = multiQuadRigidLoad(traj, params, tensions);
    
    [ref] = quadrotorLoad(traj, mQ, J, mL, l, sgnT)
    [ref] = computeLoadAttitudes(Tq, dTq, d2Tq, d3Tq, d4Tq, sgnT)
end
    
end