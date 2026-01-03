%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to compute Radar Cross Section of Dielectric Cylinder       %
% in elliptical cylinder                                              % 
% TE-Polarization                                                     %
% Written by:                                                         %
%             Seyed Mojtaba Mohammadi                                 %
%             9007704                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
E0=8.85e-12;
miu0=4*pi*1e-7;
c=3e8;
%f0=input('please insert the frequency (GHz):');
f0=1;
%phi0=input('please insert Phi0 (Radian):');
phi0=pi/4;
f0=f0*1e9;
lambda=c/f0;
%A=input('please insert A (cm):');
%A=A/100;
A=lambda/4;
%B=input('please insert B (cm):');
%B=B/100;
B=lambda/8;
f=(A^2-B^2)^0.5;
%Er1=input('please insert Er1:');
Er1=1;
%zigma1=input('please insert Zigma1:');
zigma1=0;
%miur1=input('please insert MUr1:');
miur1=1;
E1=E0*Er1-1i*zigma1/(2*pi*f0);
%Er2=input('please insert Er2:');
Er2=5;
%zigma2=input('please insert Zigma2:');
zigma2=0.1;
E2=E0*Er2-1i*zigma2/(2*pi*f0);
miu1=miu0*miur1;
%miur2=input('please insert MUr2:');
miur2=1;
miu2=miu0*miur2;
eta1=(miu1/E1)^0.5;
eta2=(miu2/E2)^0.5;
u0=acosh(A/f);
k1=2*pi*f0*(miu1*E1)^0.5;
k2=2*pi*f0*(miu2*E2)^0.5;
q1=(f^2)*(k1^2)/4;
q2=(f^2)*(k2^2)/4;
%M=input('please insert Maximum terms(M):');
M=10;
[at,as]=Amatrix(u0,phi0,q1,q2,E1,E2,M);
[bt,bs]=Bmatrix(u0,phi0,q1,q2,E1,E2,M);

Etv=[];
Etu=[];
Et1=[];
ro=10*lambda;
phi=linspace(0,2*pi,100);
phi1=phi*180/pi;
y=ro*sin(phi);
x=ro*cos(phi);
for nn=1:100
    p1=((x(nn)-f)^2+y(nn)^2)^0.5;
    p2=((x(nn)+f)^2+y(nn)^2)^0.5;
    u=acosh((p1+p2)/(2*f));
    v=acos((p2-p1)/(2*f));
    if y(nn)<0
        v=-v;
    end
    
    alfa=(1i)*(f*(((cosh(u))^2-(cos(v))^2)^0.5)*2*pi*f0*E1)^-1;
    Etv(nn)=0;
    Etu(nn)=0;

    for m=0:M
        mm=m+1;
        Etv(nn)=Etv(nn)+alfa*as(mm,1)*dMpm4(1,m,q1,u)*Spm(1,m,q1,v);
        Etu(nn)=Etu(nn)+alfa*as(mm,1)*Mpm4(1,m,q1,u)*dSpm(1,m,q1,v);
    end
    for m=1:M
        mm=m;
        Etv(nn)=Etv(nn)+alfa*bs(mm,1)*dMpm4(2,m,q1,u)*Spm(2,m,q1,v);
        Etu(nn)=Etu(nn)+alfa*bs(mm,1)*Mpm4(2,m,q1,u)*dSpm(2,m,q1,v);
    end
    
    Et1(nn)=((abs(Etv(nn)))^2+(abs(Etu(nn)))^2)^0.5;
    Et1(nn)=2*pi*ro*Et1(nn)^2/eta1^2/lambda;
end


load s2.p
plot(phi1,Et1,'b','Linewidth',3)
hold
plot(s2(:,1),s2(:,2),'r-*')
axis tight
title('Radar Cross Section of Dielectric Elliptical Cylinder, TE-Polarization','Color','b')
xlabel('\phi Angle')
ylabel('\sigma/\lambda')
legend('Series','MoM')
gtext('\phi_0=\pi/4,fre.=1[GHz], A=\lambda/4, B=\lambda/8, \epsilon_r=5., \sigma=0.1[S/m]')
print -djpeg -r300 RCS_TEz.jpg
pause(5)
close all


