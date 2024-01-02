clc
clear

%% ICRF3 catalog
%0454+844
a04 = [5*360/24 8 42.36351222]; % alpha - [deg min sec]
d04 = [84 32 04.5441733]; % delta - [deg min sec]
a04 = dms2dd(a04); d04 = dms2dd(d04); % Convert to [dec deg]
%1101−536
a1101 = [11*360/24 3 52.22168463]; 
d1101 = -1*[53 57 00.6966389];
a1101 = dms2dd(a1101); d1101 = dms2dd(d1101);
%1111+149 
a1111 = [11*360/24 13 58.69508613];
d1111 = [14 42 26.9526507];
a1111 = dms2dd(a1111); d1111 = dms2dd(d1111);
%1738+499
a17 = [17*360/24 39 27.39049431];
d17 = [49 55 03.3683385];
a17 = dms2dd(a17); d17 = dms2dd(d17);
%Spherical Coordinate of 4 Objects
obj1 = [d04,a04,1]; obj2 = [d1101,a1101,1];
obj3 = [d1111,a1111,1]; obj4 = [d17,a17,1];

%Cartesian Coordinate of 4 Objects
xi0_1 = latlonr2xyz(obj1);
xi0_2 = latlonr2xyz(obj2);
xi0_3 = latlonr2xyz(obj3);
xi0_4 = latlonr2xyz(obj4);

%Berlin
lat_berlin = [52 36 0];% deg min sec 
lon_berlin = [13 24 0];% deg min sec
lat_berlin = dms2dd(lat_berlin); lon_berlin = dms2dd(lon_berlin);
berlin = [lat_berlin,lon_berlin,1];

%polar motion
[year,month,day,mjd,xpol,ypol,dUT1,LOD,dX,dY...
s_xpol,s_ypol,s_dUT1,s_LOD,s_dX,s_dy] = textread('IAU1980_body2021.txt');

% Xpol,Ypol between 2021-06-11 & 2021-06-12 
%Row 162 - 2021-06-11
%Row 163 - 2021-06-12
xpol_inter = interp1(xpol(162:163),1:(1/24/60):2);%in min. resolution
ypol_inter = interp1(ypol(162:163),1:(1/24/60):2);%in min. resolution

% JD in 2021-06-11 & 12 with min. resolution
mm=6;yyyy = 2021;
dd=1; ut1=1:24; minute=1:60; second=0;
j = 1;
for ut1 = 1:24
    % Seperated Each Hour
    jd = gre2jd(yyyy,mm,dd,ut1,minute,second); 
    % JD of 1 min in Specific hr
    jd = jd(1:60); % Extract only Julian Date
    jd_all(j:j+59) = jd; %Insert Each Year of JD value in 1 hour
    j = j+60;
end
% The final time
lastmin = gre2jd(2021,6,12,0,0,0);
jd_all(1441) = lastmin(1);
jd_all;
% All JD in one day with min resolution


for i = 1:length(xpol_inter)

    t = (jd_all(i) - 2451545.0)/36525;
    GMST(i) = (F(jd_all(i))*86400 + 24110.54841 - 86400/2 + ...
        8640184.812866 * t + 0.093104 * t * t - ...
        6.2e-6 * t * t * t)/3600*60; % [min]
    GMST(i) = rem(GMST(i),60);

    w = pol(xpol_inter(i),ypol_inter(i));
    %Object 1
    xg_1 = ref3d(1)*rot3d(90-lat_berlin,2)*rot3d(lon_berlin,3)*...
    w*rot3d(GMST(i),3)*nut(jd_all(i))*prec(jd_all(i))*xi0_1';
    %sph_xg = xyz2latlonr(xg)
    A_1(i) = atan2(xg_1(2),xg_1(1)); %rad
    E_1(i) = atan(xg_1(3)/sqrt(xg_1(1)^2+xg_1(2)^2));%rad
    
    %Object 2
    xg_2 = ref3d(1)*rot3d(90-lat_berlin,2)*rot3d(lon_berlin,3)*...
    w*rot3d(GMST(i),3)*nut(jd_all(i))*prec(jd_all(i))*xi0_2';
    %sph_xg = xyz2latlonr(xg)
    A_2(i) = atan2(xg_2(2),xg_2(1)); %rad
    E_2(i) = atan(xg_2(3)/sqrt(xg_2(1)^2+xg_2(2)^2));%rad
    
    %Object 3
    xg_3 = ref3d(1)*rot3d(90-lat_berlin,2)*rot3d(lon_berlin,3)*...
    w*rot3d(GMST(i),3)*nut(jd_all(i))*prec(jd_all(i))*xi0_3';
    %sph_xg = xyz2latlonr(xg)
    A_3(i) = atan2(xg_3(2),xg_3(1)); %rad
    E_3(i) = atan(xg_3(3)/sqrt(xg_3(1)^2+xg_3(2)^2));%rad

    %Object 4
    xg_4 = ref3d(1)*rot3d(90-lat_berlin,2)*rot3d(lon_berlin,3)*...
    w*rot3d(GMST(i),3)*nut(jd_all(i))*prec(jd_all(i))*xi0_4';
    %sph_xg = xyz2latlonr(xg)
    A_4(i) = atan2(xg_4(2),xg_4(1)); %rad
    E_4(i) = atan(xg_4(3)/sqrt(xg_4(1)^2+xg_4(2)^2));%rad

end

A_1 = A_1*180/pi;
for m = 1:length(A_1)
    if A_1(m) < 0
        A_1(m) = A_1(m)+360;
    end
end%
A_1;
E_1 = E_1*180/pi;

A_2 = A_2*180/pi;
for n = 1:length(A_2)
    if A_2(n) < 0
        A_2(n) = A_2(n)+360;
    end
end%
A_2;
E_2 = E_2*180/pi;

A_3 = A_3*180/pi;
for o = 1:length(A_3)
    if A_3(o) < 0
        A_3(o) = A_3(o)+360;
    end
end%
A_3;
E_3 = E_3*180/pi;

A_4 = A_4*180/pi;
for q = 1:length(A_4)
    if A_4(q) < 0
        A_4(q) = A_4(q)+360;
    end
end%
A_4;
E_4 = E_4*180/pi;

% Time plot obj 1
figure(1)
xaxis = linspace(1,25,1441);
fig1 = plotyy(xaxis,A_1,xaxis,E_1);
xlabel('time (hour)')
ylabel(fig1(1),'elev (deg)')
ylabel(fig1(2),'azi (deg)')


% Time plot obj 2
figure(2)
fig2 = plotyy(xaxis,A_2,xaxis,E_2);
xlabel('time (hour)')
ylabel(fig2(1),'elev (deg)')
ylabel(fig2(2),'azi (deg)')


% Time plot obj 3
figure(3)
fig3 = plotyy(xaxis,A_3,xaxis,E_3);
xlabel('time (hour)')
ylabel(fig3(1),'elev (deg)')
ylabel(fig3(2),'azi (deg)')


% Time plot obj 4
figure(4)
fig4 = plotyy(xaxis,A_4,xaxis,E_4);
xlabel('time (hour)')
ylabel(fig4(1),'elev (deg)')
ylabel(fig4(2),'azi (deg)')

%Sky plot of 4 objs
figure(5)
all_azi = [A_1 A_2 A_3 A_4];
all_elev = [E_1 E_2 E_3 E_4];
skyplot(all_azi,all_elev,'+')
figure(6)
skyplot(all_azi,all_elev,'c+:')