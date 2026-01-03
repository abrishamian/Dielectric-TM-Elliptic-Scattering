function f4r=Mpm4(kf,m,q,x,varargin)
%     ==============================================================
%     Purpose: Compute radial Mathieu functions of the fourth kind Mcm(q,x)and Msm(q,x)
%
%     Input:   KF --- Function code
%     KF=1 for computing even function Mcm'(q,x)
%     KF=2 for computing odd function Msm'(q,x)
%     m  --- Order of Mathieu functions
%     q  ---Vector Parameter of Mathieu functions(q can be Positive or Complex number)
%     x  ---Arguments of rsdial Mathieu functions
%     Output:  f4r --- Mcm(q,x)or Msm(q,x)
%
%     ==============================================================
f1r=Mpm1(kf,m,q,x,varargin);
f2r=Mpm2(kf,m,q,x,varargin);
f4r=f1r-i.*f2r;
end