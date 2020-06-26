function [ref] = multiQuadRigidLoad(traj, params)
    

%% parameters
nQ = params.nQ;
mL = params.mL;
JL = params.JL;
rb = params.rb;
rb_hat = params.rb_hat;
g = params.g;
e1 = params.e1;
e2 = params.e2;
e3 = params.e3;

% mQ = params.mQ;
% JQ = params.JQ;
% L = params.L;


% compute the corresponding wrench and its higher derivatives in body frame
W = [mL*(traj.d2xL + g*e3); JL*traj.dOmegaL + cross(traj.OmegaL, JL*traj.OmegaL)];
dW = [mL* traj.d3xL;
      JL*traj.d2OmegaL + cross(traj.dOmegaL, JL*traj.OmegaL) ...
      + cross(traj.OmegaL, JL*traj.dOmegaL)];
  
d2W = [mL*traj.d4xL;
        JL*traj.d3OmegaL + cross(traj.d2OmegaL, JL*traj.OmegaL)...
        + 2*cross(traj.dOmegaL, JL*traj.dOmegaL) + cross(traj.OmegaL, JL*traj.d2OmegaL)];
    
d3W = [mL*traj.d5xL;
       JL*traj.d4OmegaL + cross(traj.d3OmegaL, JL*traj.OmegaL)...
       + 3*cross(traj.d2OmegaL, JL*traj.dOmegaL) + 3*cross(traj.dOmegaL, JL*traj.d2OmegaL)...
       + cross(traj.OmegaL, JL*traj.d3OmegaL)];
           
d4W = [mL*traj.d6xL;
       JL*traj.d5OmegaL + cross(traj.d4OmegaL, JL*traj.OmegaL)...
       + 4*cross(traj.d3OmegaL, JL*traj.dOmegaL) + 6*cross(traj.d2OmegaL, JL*traj.d2OmegaL)...
       + 4*cross(traj.dOmegaL, JL*traj.d3OmegaL) + cross(traj.OmegaL, JL*traj.d4OmegaL)];
% ----------------------------------------------------------------------------------------      

% << Remove the hard-code on the tension flats>>
% Tensions - Flat-Outputs (constant) % to be updated in the future
LAMBDA = [(1/4)*mL*g*vec_dot(RotY(pi/6)*e3,e1); 
                (1/4)*mL*g*vec_dot(RotY(pi/6)*e3,e2);
                (1/4)*mL*g*vec_dot(RotX(pi/6)*e3,e1); 
                (1/4)*mL*g*RotY(-pi/6)*e3];
% LAMBDA = [(1/4)*mL*g*vec_dot(RotY(0/6)*e3,e1); 
%                 (1/4)*mL*g*vec_dot(RotY(0/6)*e3,e2);
%                 (1/4)*mL*g*vec_dot(RotX(0/6)*e3,e1); 
%                 (1/4)*mL*g*RotY(0/6)*e3];
            
% << Update this to make it time varying >>
% r = [params.rb, params.rb(:,1)];
% u = cell(1,1);
% for i = 1:(params.nQ)
%    for j = 1:params.nQ
%        if i~=j
%            u{i,j} = (r(:,j)-r(:,i))/norm(r(:,j)-r(:,i));
%        end
%    end
% end
% % (FIX THIS - write code to automate following variable - data.N)
% % ----------------------------------------------------------------
% params.N = [[u{1,2};-u{1,2};zeros(3,1);zeros(3,1)],...
%             [u{1,3};zeros(3,1);-u{1,3};zeros(3,1)],...
%             [u{1,4};zeros(3,1);zeros(3,1);-u{1,4}],...
%             [zeros(3,1);u{2,3};-u{2,3};zeros(3,1)],...
%             [zeros(3,1);u{2,4};zeros(3,1);-u{2,4}],...
%             [zeros(3,1);zeros(3,1);u{3,4};-u{3,4}]];

Tq = reshape(params.Phi_pinv*W + params.N*LAMBDA, 3, nQ);            
dTq = reshape(params.Phi_pinv*dW, 3, nQ);
d2Tq = reshape(params.Phi_pinv*d2W, 3, nQ);
d3Tq = reshape(params.Phi_pinv*d3W, 3, nQ);
d4Tq = reshape(params.Phi_pinv*d4W, 3, nQ);


ref.q = zeros(3, nQ);
ref.omega = zeros(3, nQ);
ref.RQ = zeros(3,3,nQ);
ref.OmegaQ = zeros(3,nQ);
ref.dOmegaQ = zeros(3,nQ);
ref.u.f = zeros(1, nQ);
ref.u.M = zeros(3, nQ);

for i = 1:nQ
    cable = Flat2State.computeLoadAttitudes(-Tq(:,i), -dTq(:,i), -d2Tq(:,i), -d3Tq(:,i), -d4Tq(:,i), 1);
    cable.om = cross(cable.q, cable.dq);
    cable.dom = cross(cable.q, cable.d2q);

    % quadrotor position
    quad_.xQ = traj.xL + traj.RL*params.quad(i).rb - params.quad(i).l*cable.q ;
    quad_.vQ = traj.dxL + traj.dRL*params.quad(i).rb - params.quad(i).l*cable.dq ;
    quad_.aQ = traj.d2xL + traj.d2RL*params.quad(i).rb - params.quad(i).l*cable.d2q ;
    quad_.daQ = traj.d3xL + traj.d3RL*params.quad(i).rb - params.quad(i).l*cable.d3q ;
    quad_.d2aQ = traj.d4xL + traj.d4RL*params.quad(i).rb - params.quad(i).l*cable.d4q ;

    % quadrotor attitude
    F = params.quad(i).mQ*(quad_.aQ+g*e3) + Tq(:,i);
    dF = params.quad(i).mQ*(quad_.daQ) + dTq(:,i);
    d2F = params.quad(i).mQ*(quad_.d2aQ) + d2Tq(:,i);

    [r2] = Flat2State.computeQuadrotorMoment(F, dF, d2F, params.quad(i).mQ, params.quad(i).JQ);

    
    ref.q(:,i) = cable.q;
    ref.omega(:,i) = cable.om;
    ref.domega(:,i) = cable.dom;
    ref.RQ(:,:,i) = r2.R;
    ref.OmegaQ(:,i) = r2.Om;
    ref.dOmegaQ(:,i) = r2.dOm;
    ref.u.f(i) = r2.f;
    ref.u.M(:,i) = r2.M;
end

ref.xL =  traj.xL;
ref.vL =  traj.dxL;
ref.aL =  traj.d2xL;
ref.RL =  traj.RL;
ref.OmegaL = traj.OmegaL;
ref.dOmegaL = traj.dOmegaL;



end