function [P] = prec(jd)
    d = jd - 2451545.0;
    T = d/36525;
    zeta = 2306.2181 * T + 0.30188 * T^2; %sec deg.
    za = 2306.2181 * T +  1.09468 * T^2; %sec deg.
    theta = 2004.3109 * T - 0.42665 * T^2; %sec deg.
    
    % Second deg.-> Decimal deg.
    zeta = dms2dd([0 0 zeta]);
    za = dms2dd([0 0 za]);
    theta = dms2dd([0 0 theta]);

    P = rot3d(-za,3)*rot3d(theta,2)*rot3d(-zeta,3);