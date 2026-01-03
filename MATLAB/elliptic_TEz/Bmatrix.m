function [bt,bs]=Bmatrix(u0,fi0,q1,q2,miu1,miu2,M)

kf=2;

%Generate Gn & Hn
for m=1:M
    mm=m;
    bi(mm)=2*(1i^m)*Spm(kf,m,q1,fi0);
    Gn(mm,1)=bi(mm)*Mpm1(kf,m,q1,u0);
    Hn(mm,1)=(miu1^-1)*bi(mm)*dMpm1(kf,m,q1,u0);
end

%Generate Amn & Bmn & Cmn & Dmn
for n=1:M
    nn=n;
    for m=1:M
        mm=m;
        Amn(nn,mm)=Mpm1(kf,m,q2,u0)*Ipm(kf,q2,q1,m,n);
        Bmn(nn,mm)=-1*Mpm4(kf,m,q1,u0)*deltamn(m,n);
        Cmn(nn,mm)=(miu2^-1)*dMpm1(kf,m,q2,u0)*Ipm(kf,q2,q1,m,n);
        Dmn(nn,mm)=(miu1^-1)*(-1)*dMpm4(kf,m,q1,u0)*deltamn(m,n);
    end
end

A=[Amn,Bmn;Cmn,Dmn];
B=[Gn;Hn];
X=A\B;
bt=X(1:M,1);
bs=X(M+1:2*M,1);


























end