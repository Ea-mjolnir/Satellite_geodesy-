function [cat] = latlonr2xyz(sph)
    lat = sph(1); lon = sph(2); r = sph(3);
    x = r * cos(lon) * cos(lat);
    y = r * sin(lon) * cos(lat);
    z = r * sin(lat);
    [cat] = [x,y,z];
end