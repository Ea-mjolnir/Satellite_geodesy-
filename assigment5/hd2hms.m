function [hms] = hd2hms(hd)
    h = fix(hd);
    m = fix((hd-h)*60);
    s = (((hd-h)*60)-m)*60;
    [hms] = [h,m,s];
end