function d4r=dMpm4(kf,m,q,x,varargin)
%     ==============================================================
%     Purpose: Compute the derivative of radial Mathieu functions of the forrth kind Mcm'(q,x)and Msm'(q,x)
%
%     Input:   KF --- Function code
%     KF=1 for computing even function Mcm'(q,x)
%     KF=2 for computing odd function Msm'(q,x)
%     m  --- Order of Mathieu functions
%     q  ---Vector Parameter of Mathieu functions(q can be Positive or Complex number)
%     x  ---Arguments of rsdial Mathieu functions
%     Output:  f4r --- Mcm'(q,x)or Msm'(q,x)
%
%     ==============================================================
d1r=dMpm1(kf,m,q,x,varargin);
d2r=dMpm2(kf,m,q,x,varargin);
d4r=d1r-i.*d2r;
end