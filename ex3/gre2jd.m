function [jd_mjd] = gre2jd(yyyy,mm,dd,hour,minute,second)
%% gre2jd converts Gregorian calendar date and day time to Julian and Modified Julian Day
% Inputs:
%   yyyy: year
%   mm  : month
%   dd  : date
%   hour: hout
%   minute: minute
%   second: second
% Outputs:
%   [jd,mjd]: Julian and Modified Julian Day

fd = hour./24 + minute./1440 + second./86400;
my = fix((mm-14)./12);
jd = fix((1461.*(yyyy + 4800 + my))./4) + ...
    fix((367.*((mm-2)-(12.*my)))./12) - ...
    fix((3.*((yyyy + 4900 +my)./100))./4) + ...
    dd - 32075.5 + fd;
mjd = jd-2400000.5;
jd_mjd = [jd,mjd];