function [W] = pol(xp,yp) 
    
    %xp,yp [sec. deg.] -> [decimal. deg.]
    xp = dms2dd([0 0 xp]); %Decimal deg.
    yp = dms2dd([0 0 yp]); %Decimal deg.
    
    W = rot3d(-xp,2)*rot3d(-yp,1);
