function [ref] = quadrotorLoad(traj, mQ, J, mL, l, sgnT)
% ------------------------------------------------------------
% function [ref] = quadrotorLoad(traj, mQ, J, mL, l)
% Diffferential flatness for quadrotorLoad in 3D
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
%         xL, vL, aL, q, dq, d2q, om, dom, f
%         xQ, vQ, aQ, R, dR, d2R, Om, dOm, M
% ------------------------------------------------------------

g = 9.81;
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

xL = traj.x;
vL = traj.dx;
aL = traj.d2x;
daL = traj.d3x;
d2aL = traj.d4x;
d3aL = traj.d5x;
d4aL = traj.d6x;


% load attitude
Tq = -mL*(aL + g*e3) ;
dTq = -mL*daL;
d2Tq = -mL*d2aL;
d3Tq = -mL*d3aL;
d4Tq = -mL*d4aL;
[ref] = Flat2State.computeLoadAttitudes(Tq, dTq, d2Tq, d3Tq, d4Tq, sgnT);    
ref.om = cross(ref.q, ref.dq);
    
% quadrotor position
ref.xQ = xL - l*ref.q ;
ref.vQ = vL - l*ref.dq ;
ref.aQ = aL - l*ref.d2q ;
ref.daQ = daL - l*ref.d3q ;
ref.d2aQ = d2aL - l*ref.d4q ;

% quadrotor attitude
F = mQ*(ref.aQ+g*e3) - Tq;
dF = mQ*(ref.daQ) -dTq;
d2F = mQ*(ref.d2aQ) -d2Tq;

[r2] = Flat2State.computeQuadrotorMoment(F, dF, d2F, mQ, J);

ref.xL =  xL;
ref.vL =  vL;
ref.aL =  aL;
ref.R = r2.R;
ref.Om = r2.Om;
ref.dOm = r2.dOm;
ref.M = r2.M;
ref.f = r2.f;
end