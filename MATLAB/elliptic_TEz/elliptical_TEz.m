%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to compute distribution of electric field along X & Y axis, %
% in elliptical cylinder                                              % 
% TE-Polarization                                                     %
% Written by:                                                         %
%             Seyed Mojtaba Mohammadi                                 %
%             9007704                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
E0=8.85e-12;
miu0=4*pi*1e-7;
c=3e8;
% %fi0=input('please insert Phi0 (Radian):');
fi0=pi/4;

f0=1;
f0=f0*1e9;
lambda=c/f0;
% A=input('please insert A (cm):');
A=lambda/4;
% B=input('please insert B (cm):');
B=lambda/8;
f=(A^2-B^2)^0.5;
%Er1=input('please insert Er1:');
Er1=1;
%zigma1=input('please insert Zigma1:');
zigma1=0.0;
%miur1=input('please insert MUr1:');
miur1=1;
E1=E0*Er1-1i*zigma1/(2*pi*f0);
%Er2=input('please insert Er2:');
Er2=5.;
%zigma2=input('please insert Zigma2:');
zigma2=0.1;
E2=E0*Er2-1i*zigma2/(2*pi*f0);
miu1=miu0*miur1;
%miur2=input('please insert MUr2:');
miur2=1;
miu2=miu0*miur2;
u0=acosh(A/f);
k1=2*pi*f0*(miu1*E1)^0.5;
k2=2*pi*f0*(miu2*E2)^0.5;
q1=(f^2)*(k1^2)/4;
q2=(f^2)*(k2^2)/4;
%  M=input('please insert Maximum terms(M):');
M=8;
[at,as]=Amatrix(u0,fi0,q1,q2,E1,E2,M);
[bt,bs]=Bmatrix(u0,fi0,q1,q2,E1,E2,M);

% X-axis
Etv=[];
Etu=[];
Et1=[];
y=0;
x=linspace(-A,A,100);
xx=x*100;
for nn=1:100
    p1=((x(nn)-f)^2+y^2)^0.5;
    p2=((x(nn)+f)^2+y^2)^0.5;
    u=acosh((p1+p2)/(2*f));
    v=acos((p2-p1)/(2*f));
    if y<0
        v=-v;
    end
    
    alfa=(1i)*(f*(((cosh(u))^2-(cos(v))^2)^0.5)*2*pi*f0*E2)^-1;
    Etv(nn)=0;
    Etu(nn)=0;

    for m=0:M
        mm=m+1;
        Etv(nn)=Etv(nn)+alfa*at(mm,1)*dMpm1(1,m,q2,u)*Spm(1,m,q2,v);
        Etu(nn)=Etu(nn)+alfa*at(mm,1)*Mpm1(1,m,q2,u)*dSpm(1,m,q2,v);
    end
    for m=1:M
        mm=m;
        Etv(nn)=Etv(nn)+alfa*bt(mm,1)*dMpm1(2,m,q2,u)*Spm(2,m,q2,v);
        Etu(nn)=Etu(nn)+alfa*bt(mm,1)*Mpm1(2,m,q2,u)*dSpm(2,m,q2,v);
    end
    
    Et1(nn)=((abs(Etv(nn)))^2+(abs(Etu(nn)))^2)^0.5;
end



 % Y-axis
Etv=[];
Etu=[];
Et2=[];
x=0;
y=linspace(-B,B,100);
yy=y*100;
for nn=1:100
    p1=((x-f)^2+y(nn)^2)^0.5;
    p2=((x+f)^2+y(nn)^2)^0.5;
    u=acosh((p1+p2)/(2*f));
    v=acos((p2-p1)/(2*f));
    if y(nn)<0
        v=-v;
    end
    
    alfa=(1i)*(f*(((cosh(u))^2-(cos(v))^2)^0.5)*2*pi*f0*E2)^-1;
    Etv(nn)=0;
    Etu(nn)=0;
    for m=0:M
        mm=m+1;
        Etv(nn)=Etv(nn)+alfa*at(mm,1)*dMpm1(1,m,q2,u)*Spm(1,m,q2,v);
        Etu(nn)=Etu(nn)+alfa*at(mm,1)*Mpm1(1,m,q2,u)*dSpm(1,m,q2,v);
    end
    for m=1:M
        mm=m;
        Etv(nn)=Etv(nn)+alfa*bt(mm,1)*dMpm1(2,m,q2,u)*Spm(2,m,q2,v);
        Etu(nn)=Etu(nn)+alfa*bt(mm,1)*Mpm1(2,m,q2,u)*dSpm(2,m,q2,v);
    end
    
    Et2(nn)=((abs(Etv(nn)))^2+(abs(Etu(nn)))^2)^0.5;
    
end
%

plot(xx,Et1,'b-',yy,Et2,'r','Linewidth',3)
hold
load sx.p
load sy.p
plot(sx(:,1),sx(:,2),'r*');
plot(sy(:,1),sy(:,2),'b*');
legend('Series : x','Series : y','MoM : x ','MoM : y',2)
title('Distribution of electric filed along X & Y axis, TE-Polarization')
xlabel('Distance [cm]')
ylabel('Electric Field [V/m]')
gtext('\phi_0=\pi/4,fre.=1[GHz], a=\lambda/2,b=\lambda/4 \epsilon_r=5., \sigma=0.1[S/m]')
pause(7)
print -djpeg -r300 elte5.jpg
close all
