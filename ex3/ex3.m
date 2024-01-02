clc
clear

clc
clear
% vernaleqinox
alpha = 0; delta = 0;
vernaleqinox = [delta,alpha,1];
xi0 = latlonr2xyz(vernaleqinox);  %(rad,rad,m)-->(m,m,m)

[year,month,day,mjd,xpol,ypol,dUT1,LOD,dX,dY...
s_xpol,s_ypol,s_dUT1,s_LOD,s_dX,s_dy] = textread('IAU1980_body2021.txt');

% Xpol,Ypol between 2021-06-11 & 2021-06-12
%Row 162 - 2021-06-11
%Row 163 - 2021-06-12
xpol_inter = interp1(1:length(xpol),xpol,xpol(162:163),'linear');
ypol_inter = interp1(1:length(ypol),ypol,ypol(162:163),'linear');

%W(Xp,Yp)
for i = 1:length(xpol)
    w = pol(xpol(i),ypol(i));
    xiw = w*xi0';
    % (m,m,m) --> (rad,rad,m)
    sph_iw = xyz2latlonr(xiw); 
    % [delta,alpha,1] - [rad,rad,m]
    
    delta2(i) = sph_iw(1)*180/pi; % degree
    alpha2(i) = sph_iw(2)*12/pi; % hour angle
end
delta2;
alpha2;