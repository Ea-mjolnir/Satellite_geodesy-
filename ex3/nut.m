function [N] = nut(jd)
    d = jd - 2451545.0;
    T = d/36525;

    f1 = 125 - 0.05295 * d; %Decimal deg.
    f2 = 200.9 + 1.97129 * d; %Decimal deg.
    ea = 84381.448 - 46.8150*T; %Sec deg.
    dpsi = -0.0048*sind(f1) - 0.0004*sind(f2); %Decimal deg.
    de = 0.0026*cosd(f1) - 0.0002*cosd(f2); %Decimal deg.
    
    % Second deg.-> Decimal deg.
    ea = dms2dd([0 0 ea]);

    N = rot3d(-ea-de,1)*rot3d(-dpsi,3)*rot3d(ea,1);
