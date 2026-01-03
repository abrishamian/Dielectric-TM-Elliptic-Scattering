%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CORRELATION FACTOR    
%   y = Ipm(kf,q1,q2,m,n)   [p,m = e,o (even,odd)]
%
%   INPUTS: 
%     KF=1 for  Icc(q,xr)
%     KF=2 for  Iss(q,xr)  
%     m,n   --- Order of Mathieu functions
%     q1,q2   --- Vector Parameter of Mathieu functions(q can be Positive or Complex number)
%   Ipm FUNCTION CALL
function y = Ipm(kf,q1,q2,m,n)


tt=0;
if(kf == 1&&m == 2*fix(m./2))
    KF=1;
    if n ~= 2*fix(n./2)
        y=0;
        tt=1;
    end
end;
if(kf == 1&&m ~= 2*fix(m./2))
    KF=2;
    if n == 2*fix(n./2)
        y=0;
        tt=1;
    end
end;
if(kf == 2&&m ~= 2*fix(m./2))
    KF=3;
    if n == 2*fix(n./2)
        y=0;
        tt=1;
    end
end;
if(kf == 2&&m == 2*fix(m./2))
    KF=4;
    if n ~= 2*fix(n./2)
        y=0;
        tt=1;
    end
end;
if tt==0
    
Apm=fcoef1(KF,m,q1);
AApm=fcoef1(KF,n,q2);
ncoeffs=min(length(Apm),length(AApm));

if KF == 1  

    A0=Apm(1);     AA0=AApm(1); 
    y = 2*A0*AA0;
    for j = 2:ncoeffs
        y = y +  Apm(j) * AApm(j);
    end
            
    
else
    y=0;      
    for j = 1:ncoeffs
        y = y +  Apm(j) * AApm(j);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
    