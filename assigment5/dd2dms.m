function [dms] = dd2dms(dd)
    d = fix(dd);
    m = fix((dd-d)*60);
    s = (((dd-d)*60)-m)*60;
    [dms] = [d,m,s];
end