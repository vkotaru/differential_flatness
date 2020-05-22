function [ref] = quadrotor(traj, mQ, J)
% ------------------------------------------------------------
% function [ref] = quadrotor(traj, mQ, J)
% Diffferential flatness for quadrotor in 3D
% 
% Inputs: traj, mQ, J
%         traj: is struct with fields 'x', 'dx', 'd2x',.., Psi, dPsi, d2Psi
%         containing the flat-outputs and their derivatives
%         
%         mQ: mass of the quadrotor
%          J: inertia of the quadrotor in the body-frame
%
% Output: ref
%         ref: is struct with fields
%         xQ, vQ, aQ, R, dR, d2R, Om, dOm, M
% ------------------------------------------------------------

g = 9.81;
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

xQ = traj.x;
vQ = traj.dx;
aQ = traj.d2x;
daQ = traj.d3x;
d2aQ = traj.d4x;


F = mQ*(aQ+g*e3);
dF = mQ*(daQ);
d2F = mQ*(d2aQ);

[ref] = Flat2State.computeQuadrotorMoment(F, dF, d2F, mQ, J);
ref.xQ = xQ;
ref.vQ = vQ;
ref.aQ = aQ;

end