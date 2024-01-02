function dd = dms2dd(dms)
    if dms(1) >= 0
       if sign(dms(2)) == -1
       dd = (-1)*(abs(dms(1))+abs(dms(2)/60)+abs(dms(3)/3600));
       elseif sign(dms(3)) == -1
       dd = (-1)*(abs(dms(1))+abs(dms(2)/60)+abs(dms(3)/3600));
       else
       dd = dms(1) + dms(2)/60 + dms(3)/3600;
       end
    else 
    dd = (-1)*(abs(dms(1))+abs(dms(2)/60)+abs(dms(3)/3600));
end