function [ref] = computeQuadrotorMoment(F, dF, d2F, varargin)
    
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

g = 9.81;

% << TODO remove this hard-code in the future >> 
b1d = e1 ;
db1d = zeros(3,1) ;
d2b1d = zeros(3,1) ;

mQ = varargin{1};
J = varargin{2};

%% R
	
% F = mQ*(aQ+g*e3);
norm_fb3 = norm(F);
f = norm_fb3 ;
b3 = F / norm_fb3 ;
b3_b1d = vec_cross(b3, b1d) ;
norm_b3_b1d = norm(b3_b1d) ;
b1 = - vec_cross(b3, b3_b1d) / norm_b3_b1d ;
b2 = vec_cross(b3, b1) ;
R = [b1 b2 b3] ;

%% dR

% dF = mQ*(daQ);
dnorm_fb3 = vec_dot(F, dF) / norm_fb3 ;
db3 = (dF*norm_fb3 - F*dnorm_fb3) / norm_fb3^2 ;
db3_b1d = vec_cross(db3, b1d) + vec_cross(b3, db1d) ;
dnorm_b3_b1d = vec_dot(b3_b1d, db3_b1d) / norm_b3_b1d ;
db1 = (-vec_cross(db3,b3_b1d)-vec_cross(b3,db3_b1d) - b1*dnorm_b3_b1d) / norm_b3_b1d ;
db2 = vec_cross(db3, b1) + vec_cross(b3, db1) ;
dR = [db1 db2 db3] ;
Omega = vee_map(R'*dR) ;

%% d2R

% d2F = mQ*(d2aQ);
d2norm_fb3 = (vec_dot(dF, dF)+vec_dot(F, d2F) - dnorm_fb3*dnorm_fb3) / norm_fb3 ;
d2b3 = ( (d2F*norm_fb3+dF*dnorm_fb3 - dF*dnorm_fb3-F*d2norm_fb3)*norm_fb3^2 - db3*norm_fb3^2*2*norm_fb3*dnorm_fb3 ) / norm_fb3^4 ;
d2b3_b1d = vec_cross(d2b3, b1d)+vec_cross(db3, db1d) + vec_cross(db3, db1d)+vec_cross(b3, d2b1d) ;
d2norm_b3_b1d = ( (vec_dot(db3_b1d,db3_b1d)+vec_dot(b3_b1d,d2b3_b1d))*norm_b3_b1d - vec_dot(b3_b1d, db3_b1d)*dnorm_b3_b1d ) / norm_b3_b1d^2 ;
d2b1 = ( (-vec_cross(d2b3,b3_b1d)-vec_cross(db3,db3_b1d) - vec_cross(db3,db3_b1d)-vec_cross(b3,d2b3_b1d) - db1*dnorm_b3_b1d-b1*d2norm_b3_b1d )*norm_b3_b1d - db1*norm_b3_b1d*dnorm_b3_b1d ) / norm_b3_b1d^2 ;
d2b2 = vec_cross(d2b3, b1)+vec_cross(db3, db1) + vec_cross(db3, db1)+vec_cross(b3, d2b1) ;
d2R = [d2b1 d2b2 d2b3] ;
dOmega = vee_map( dR'*dR + R'*d2R ) ; %vee_map( dR'*dR + R'*d2R, true ) ;

M = J*dOmega + vec_cross(Omega, J*Omega) ;

%% FINAL DESIRED TRAJECTORY
% ========================
% Quad Attitude
ref.R = R;
ref.Om = Omega;
ref.dOm = dOmega;
ref.dR = dR;
ref.d2R = d2R;

% inputs
ref.f = f;
ref.M = M;

end
