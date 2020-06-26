function [ref] = computeLoadAttitudes(Tq, dTq, d2Tq, d3Tq, d4Tq, sgnT)
% default = -1
    
norm_Tq = norm(Tq);
q = sgnT*Tq/norm_Tq ;
    
dnorm_Tq = 1/norm_Tq * vec_dot(Tq, dTq) ;
dq = (sgnT*dTq - q*dnorm_Tq) / norm_Tq ;

    
d2norm_Tq = ( vec_dot(dTq, dTq) + vec_dot(Tq, d2Tq) - dnorm_Tq^2 ) / norm_Tq ;
d2q = (sgnT*d2Tq - dq*dnorm_Tq - q*d2norm_Tq - dq*dnorm_Tq) / norm_Tq ;
    
d3norm_Tq = ( 2*vec_dot(d2Tq, dTq) + vec_dot(dTq, d2Tq)+vec_dot(Tq, d3Tq) - 3*dnorm_Tq*d2norm_Tq) / norm_Tq ;
d3q = (sgnT*d3Tq - d2q*dnorm_Tq-dq*d2norm_Tq - dq*d2norm_Tq-q*d3norm_Tq - d2q*dnorm_Tq-dq*d2norm_Tq - d2q*dnorm_Tq) / norm_Tq ;
    
d4norm_Tq = ( 2*vec_dot(d3Tq, dTq)+2*vec_dot(d2Tq, d2Tq) + vec_dot(d2Tq, d2Tq)+vec_dot(dTq, d3Tq) + vec_dot(dTq, d3Tq)+vec_dot(Tq, d4Tq) - 3*d2norm_Tq^2-3*dnorm_Tq*d3norm_Tq ...
        - d3norm_Tq*dnorm_Tq) / norm_Tq ;
d4q = (sgnT*d4Tq - d3q*dnorm_Tq-d2q*d2norm_Tq - d2q*d2norm_Tq-dq*d3norm_Tq - d2q*d2norm_Tq-dq*d3norm_Tq - dq*d3norm_Tq-q*d4norm_Tq ...
        - d3q*dnorm_Tq-d2q*d2norm_Tq - d2q*d2norm_Tq-dq*d3norm_Tq - d3q*dnorm_Tq-d2q*d2norm_Tq - d3q*dnorm_Tq ) / norm_Tq ;
  

ref.q = q;
ref.dq = dq;
ref.d2q = d2q;
ref.d3q = d3q;
ref.d4q = d4q;
    
end