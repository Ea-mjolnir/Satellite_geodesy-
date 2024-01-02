function [sph] = xyz2latlonr(cat)
    x = cat(1); y = cat(2); z = cat(3);
    lat = atan(z/sqrt(x^2+y^2));
    lon = atan2(y,x);
    r = sqrt(x^2+y^2+z^2);
    [sph] = [lat,lon,r];
end