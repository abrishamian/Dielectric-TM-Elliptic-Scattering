close;
clear all

load cx.p
load cy.p
plot(cx(:,1),cx(:,2),'r','Linewidth',3);
hold
plot(cy(:,1),cy(:,2),'b','Linewidth',3);

load sx.p
load sy.p
plot(sx(:,1),sx(:,2),'b-*');
plot(sy(:,1),sy(:,2),'r-*');
legend('Series : x','Series : y','MoM : x ','MoM : y',3)
title('Distribution of electric filed along X & Y axis, TE-Polarization')
xlabel('Distance [cm]')
ylabel('Electric Field [V/m]')
 gtext('f=1[GHz], a=\lambda/4,b=\lambda/4 \epsilon_r=5., \sigma=0.1[S/m]')
pause(4)
print -djpeg -r300 elte.jpg
close all