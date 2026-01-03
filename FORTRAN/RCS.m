       % matlab file for your plot ( remove c and use it in matlab)

       clear all
       close
       load s1.p
       load s2.p
       plot(s1(:,1),s1(:,2),'linewidth',2.5);
       hold;
       plot(s2(:,1),s2(:,2),'r*');
       xlabel('\Phi Angle')
       ylabel(' \sigma/\lambda')
       title('Radar Cross Section of Lossy Dielectric Elliptic Cylinder,TM-Polarization')
       axis tight
       legend('Series','MoM')
       gtext('\phi_0=45^o;    f=1[GHz]')
       gtext('A=\lambda/2;    B=\lambda/4 ')
       gtext('\epsilon_{r}=5.;   \sigma= 0.1[S/m]')
       pause(5)
    %   print -deps2 -r300 RCS.EPS
       print -djpeg -r150 RCS.jpg
       close all

