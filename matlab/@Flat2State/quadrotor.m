function [ref] = quadrotor(traj, mQ, J)
% Diffferential flatness for quadrotor in 3D
% 
% Inputs: traj, mQ, J
%         traj: is struct with fields 'x', 'dx', 'd2x',.., Psi, dPsi, d2Psi
%         containing the flat-outputs and their derivatives
%         
%         mQ: mass of the quadrotor
%          J: inertia of the quadrotor in the body-frame

g = 9.81;

xQ = traj.x;
vQ = traj.dx;
aQ = traj.d2x;
daQ = traj.d3x;
d2aQ = traj.d4x;

[ref] = Flat2State.computeQuadrotorMoment(aQ, daQ, d2aQ, mQ, J);
ref.xQ = xQ;
ref.vQ = vQ;
ref.aQ = aQ;

end