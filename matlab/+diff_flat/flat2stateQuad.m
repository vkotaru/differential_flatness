function [ref] = flat2stateQuad(obj,flats)
%% 
% Function to calculate the feedforward reference states and inputs for a
% given position trajectory using differential flatness



mQ = obj.m;
J = obj.J;
g = obj.g;

e3 = obj.e3;
b1d = obj.e1;
db1d = zeros(3,1);
d2b1d = zeros(3,1);

axQ = flats.d2x;
daxQ = flats.d3x;
d2axQ = flats.d4x;

fb3 = mQ*(axQ+g*e3);
norm_fb3 = norm(fb3) ;
f = norm_fb3 ;
b3 = fb3 / norm_fb3 ;
b3_b1d = cross(b3, b1d) ;
norm_b3_b1d = norm(b3_b1d) ;
b1 = - cross(b3, b3_b1d) / norm_b3_b1d ;
b2 = cross(b3, b1) ;
R = [b1 b2 b3] ;

dfb3 = mQ*(daxQ);
dnorm_fb3 = dot(fb3, dfb3) / norm_fb3 ;
db3 = (dfb3*norm_fb3 - fb3*dnorm_fb3) / norm_fb3^2 ;
db3_b1d = cross(db3, b1d) + cross(b3, db1d) ;
dnorm_b3_b1d = dot(b3_b1d, db3_b1d) / norm_b3_b1d ;
db1 = (-cross(db3,b3_b1d)-cross(b3,db3_b1d) - b1*dnorm_b3_b1d) / norm_b3_b1d ;
db2 = cross(db3, b1) + cross(b3, db1) ;
dR = [db1 db2 db3] ;
Omega = vee(R'*dR) ;

d2fb3 = mQ*(d2axQ);
d2norm_fb3 = (dot(dfb3, dfb3)+dot(fb3, d2fb3) - dnorm_fb3*dnorm_fb3) / norm_fb3 ;
d2b3 = ( (d2fb3*norm_fb3+dfb3*dnorm_fb3 - dfb3*dnorm_fb3-fb3*d2norm_fb3)*norm_fb3^2 - db3*norm_fb3^2*2*norm_fb3*dnorm_fb3 ) / norm_fb3^4 ;
d2b3_b1d = cross(d2b3, b1d)+cross(db3, db1d) + cross(db3, db1d)+cross(b3, d2b1d) ;
d2norm_b3_b1d = ( (dot(db3_b1d,db3_b1d)+dot(b3_b1d,d2b3_b1d))*norm_b3_b1d - dot(b3_b1d, db3_b1d)*dnorm_b3_b1d ) / norm_b3_b1d^2 ;
d2b1 = ( (-cross(d2b3,b3_b1d)-cross(db3,db3_b1d) - cross(db3,db3_b1d)-cross(b3,d2b3_b1d) - db1*dnorm_b3_b1d-b1*d2norm_b3_b1d )*norm_b3_b1d - db1*norm_b3_b1d*dnorm_b3_b1d ) / norm_b3_b1d^2 ;
d2b2 = cross(d2b3, b1)+cross(db3, db1) + cross(db3, db1)+cross(b3, d2b1) ;
d2R = [d2b1 d2b2 d2b3] ;
dOmega = vee( dR'*dR + R'*d2R ) ; %vee_map( dR'*dR + R'*d2R, true ) ;
    
M = J*dOmega + cross(Omega, J*Omega) ;
F = fb3;


ref.x = flats.x;
ref.v = flats.dx;
ref.a = flats.d2x;
ref.da = flats.d3x;
ref.d2a = flats.d4x;

ref.R = R;
ref.Omega = Omega;
ref.dOmega = dOmega;
ref.F =F;
ref.M = M;
ref.f = dot(F,R*obj.e3);

end


