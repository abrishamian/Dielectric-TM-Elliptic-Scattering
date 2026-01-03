close all;
clear all
clc
load cx.p
load cy.p
load sx.p
load sy.p

plot(cx(:,1),cx(:,2),'linewidth',2.5);
hold;
plot(sx(:,1),sx(:,2),'r*');
%
plot(cy(:,1),cy(:,2),'r-.','linewidth',2.5);
plot(sy(:,1),sy(:,2),'v');
legend('Series : x','Series : y','MoM : x ','MoM : y',2)
title('Distribution of electric filed along X & Y axis, TM-Polarization')
xlabel('Distance [cm]')
ylabel('Electric Field [V/m]')
axis tight
gtext('\phi_0=45^o;   f=1[GHz]')
gtext('A=\lambda/2;  B=\lambda/4 ')
gtext('\epsilon_{r}= 5.   \sigma=0.1[S/m]')
print -djpeg -r150 ELXY.jpg
pause(5)
close all