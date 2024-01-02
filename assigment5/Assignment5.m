clc
clear
format long 
%%% 6 Satellites %%%
 %  i inclination (deg), Right Ascension of the Ascending Node (deg) 
 %  e Eccentricity, Argument of Perigee (deg), M Mean Anomaly (deg) 
 %  n Mean Motion (revolutions per day), t (date, time UT1)
ISS = [51.7, 19.1, 0.00, 17.6, 10.4, 15.50];
t_ISS = datetime(2022,01,14,04,49,0);

Sentinel1B = [98.2, 23.5, 0.00, 80.6, 279.5, 14.59];
t_Sentinel1B = datetime(2022,01,13,22,15,0);

Molniya1T = [63.6, 198.3, 0.62, 299.8, 13.0, 3.19]; 
t_Molniya1T = datetime(2022,01,12,18,29,0);

Galileo9 = [55.7, 266.4, 0.00, 14.0, 346.1, 1.71];
t_Galileo9 = datetime(2022,01,09,02,07,0);

Beidou3 = [60.1, 55.8, 0.00, 192.4, 355.6, 1.00]; 
t_Beidou3 = datetime(2022,01,13,01,41,0);

Meteosat11 = [0.3, 325.8, 0.00, 117.3, 338.5, 1.00]; 
t_Meteosat11 = datetime(2022,01,13,20,35,0);

%% 1
%%% Find Tsat (Tsat = 1 / nsat) %%%
T_ISS_d = 1/ISS(6); %[d]
T_ISS_h = T_ISS_d*24; %[h in decimal degree]
T_ISS_hms = hd2hms(T_ISS_h); %[h in h min sec)]
T_ISS_hm = T_ISS_hms(1:2); %Round to min

T_Sentinel1B_d = 1/Sentinel1B(6); %[d]
T_Sentinel1B_h = T_Sentinel1B_d*24; %[h in decimal degree]
T_Sentinel1B_hms = hd2hms(T_Sentinel1B_h); %[h in h min sec)]
T_Sentinel1B_hm = T_Sentinel1B_hms(1:2); %Round to min

T_Molniya1T_d = 1/Molniya1T(6); %[d]
T_Molniya1T_h = T_Molniya1T_d*24; %[h in decimal degree]
T_Molniya1T_hms = hd2hms(T_Molniya1T_h); %[h in h min sec)]
T_Molniya1T_hm = T_Molniya1T_hms(1:2); %Round to min

T_Galileo9_d = 1/Galileo9(6); %[d]
T_Galileo9_h = T_Galileo9_d*24; %[h in decimal degree]
T_Galileo9_hms = hd2hms(T_Galileo9_h); %[h in h min sec)]
T_Galileo9_hm = T_Galileo9_hms(1:2); %Round to min

T_Beidou3_d = 1/Beidou3(6); %[d]
T_Beidou3_h = T_Beidou3_d*24; %[h in decimal degree]
T_Beidou3_hms = hd2hms(T_Beidou3_h); %[h in h min sec)]
T_Beidou3_hm = T_Beidou3_hms(1:2); %Round to min


T_Meteosat11_d = 1/Meteosat11(6); %[d]
T_Meteosat11_h = T_Meteosat11_d*24; %[h in decimal degree]
T_Meteosat11_hms = hd2hms(T_Meteosat11_h); %[h in h min sec)]
T_Meteosat11_hm = T_Meteosat11_hms(1:2); %Round to min

%% 2
%%% Find major axis a %%%
 % gravitational constant
GM = 398600.44; %[km^3/s^2]
 % a = ( (Tsat^2 * GM) / (4*(pi^2)) )^1/3 

T_ISS_sec = T_ISS_h*3600; %[sec in decimal degree]
a_ISS_km = ((T_ISS_sec^2 * GM) / (4*(pi^2)))^(1/3); %[km]
a_ISS_m = round(a_ISS_km*1000); %[m]
 % Perigee & Apogee
PerigeeHeight_ISS = a_ISS_m*(1-ISS(3)); %[m]
ApogeeHeight_ISS = a_ISS_m*(1+ISS(3));  %[m]

T_Sentinel1B_sec = T_Sentinel1B_h*3600; %[sec in decimal degree]
a_Sentinel1B_km = ((T_Sentinel1B_sec^2 * GM) / (4*(pi^2)))^(1/3); %[km]
a_Sentinel1B_m = round(a_Sentinel1B_km*1000); %[m]
 % Perigee & Apogee
PerigeeHeight_Sentinel1B = a_Sentinel1B_m*(1-Sentinel1B(3)); %[m]
ApogeeHeight_Sentinel1B = a_Sentinel1B_m*(1+Sentinel1B(3));  %[m]

T_Molniya1T_sec = T_Molniya1T_h*3600; %[sec in decimal degree]
a_Molniya1T_km = ((T_Molniya1T_sec^2 * GM) / (4*(pi^2)))^(1/3); %[km]
a_Molniya1T_m = round(a_Molniya1T_km*1000); %[m]
 % Perigee & Apogee
PerigeeHeight_Molniya1T = a_Molniya1T_m*(1-Molniya1T(3)); %[m]
ApogeeHeight_Molniya1T = a_Molniya1T_m*(1+Molniya1T(3));  %[m]

T_Galileo9_sec = T_Galileo9_h*3600; %[sec in decimal degree]
a_Galileo9_km = ((T_Galileo9_sec^2 * GM) / (4*(pi^2)))^(1/3); %[km]
a_Galileo9_m = round(a_Galileo9_km*1000); %[m]
 % Perigee & Apogee
PerigeeHeight_Galileo9 = a_Galileo9_m*(1-Galileo9(3)); %[m]
ApogeeHeight_Galileo9 = a_Galileo9_m*(1+Galileo9(3));  %[m]

T_Beidou3_sec = T_Beidou3_h*3600; %[sec in decimal degree]
a_Beidou3_km = ((T_Beidou3_sec^2 * GM) / (4*(pi^2)))^(1/3); %[km]
a_Beidou3_m = round(a_Beidou3_km*1000); %[m]
 % Perigee & Apogee
PerigeeHeight_Beidou3 = a_Beidou3_m*(1-Beidou3(3)); %[m]
ApogeeHeight_Beidou3 = a_Beidou3_m*(1+Beidou3(3));  %[m]

T_Meteosat11_sec = T_Meteosat11_h*3600; %[sec in decimal degree]
a_Meteosat11_km = ((T_Meteosat11_sec^2 * GM) / (4*(pi^2)))^(1/3); %[km]
a_Meteosat11_m = round(a_Meteosat11_km*1000); %[m]
 % Perigee & Apogee
PerigeeHeight_Meteosat11 = a_Meteosat11_m*(1-Meteosat11(3)); %[m]
ApogeeHeight_Meteosat11 = a_Meteosat11_m*(1+Meteosat11(3));  %[m]

%% 3
i = 1:1440;
t_min(i) = (i-1); %min
t_d(i) = t_min(i)/1440; %d or day

Mi_ISS(i) = t_d(i) * ISS(6) * (2*pi); %rad
Msati_ISS(i) = ISS(5)/180*pi + Mi_ISS(i); %rad
Esati_ISS = kepler_e_9(Msati_ISS,ISS(3)); %rad
rsati_ISS_km = a_ISS_km*(1-ISS(3)*cos(Esati_ISS)); %km
vsati_ISS = atan2((sqrt(1-(ISS(3)^2)))*sin(Esati_ISS),(cos(Esati_ISS)-ISS(3))); %rad
rorbsati_ISS = rsati_ISS_km .* [cos(vsati_ISS);sin(vsati_ISS);zeros(1,1440)];%km
rECI_ISS = rot3d(-ISS(2),3)*rot3d(-ISS(1),1)*rot3d(-ISS(4),3)*rorbsati_ISS; %km

Mi_Sentinel1B(i) = t_d(i) * Sentinel1B(6) * (2*pi); %rad
Msati_Sentinel1B(i) = Sentinel1B(5)/180*pi + Mi_Sentinel1B(i); %rad
Esati_Sentinel1B = kepler_e_9(Msati_Sentinel1B,Sentinel1B(3)); %rad
rsati_Sentinel1B_km = a_Sentinel1B_km*(1-Sentinel1B(3)*cos(Esati_Sentinel1B)); %km
vsati_Sentinel1B = atan2((sqrt(1-(Sentinel1B(3)^2)))*sin(Esati_Sentinel1B),(cos(Esati_Sentinel1B)-Sentinel1B(3))); %rad
rorbsati_Sentinel1B = rsati_Sentinel1B_km .* [cos(vsati_Sentinel1B);sin(vsati_Sentinel1B);zeros(1,1440)];%km
rECI_Sentinel1B = rot3d(-Sentinel1B(2),3)*rot3d(-Sentinel1B(1),1)*rot3d(-Sentinel1B(4),3)*rorbsati_Sentinel1B; %km

Mi_Molniya1T(i) = t_d(i) * Molniya1T(6) * (2*pi); %rad
Msati_Molniya1T(i) = Molniya1T(5)/180*pi + Mi_Molniya1T(i); %rad
Esati_Molniya1T = kepler_e_9(Msati_Molniya1T,Molniya1T(3)); %rad
rsati_Molniya1T_km = a_Molniya1T_km*(1-Molniya1T(3)*cos(Esati_Molniya1T)); %km
vsati_Molniya1T = atan2((sqrt(1-(Molniya1T(3)^2)))*sin(Esati_Molniya1T),(cos(Esati_Molniya1T)-Molniya1T(3))); %rad
rorbsati_Molniya1T = rsati_Molniya1T_km .* [cos(vsati_Molniya1T);sin(vsati_Molniya1T);zeros(1,1440)];%km
rECI_Molniya1T = rot3d(-Molniya1T(2),3)*rot3d(-Molniya1T(1),1)*rot3d(-Molniya1T(4),3)*rorbsati_Molniya1T; %km

Mi_Galileo9(i) = t_d(i) * Galileo9(6) * (2*pi); %rad
Msati_Galileo9(i) = Galileo9(5)/180*pi + Mi_Galileo9(i); %rad
Esati_Galileo9 = kepler_e_9(Msati_Galileo9,Galileo9(3)); %rad
rsati_Galileo9_km = a_Galileo9_km*(1-Galileo9(3)*cos(Esati_Galileo9)); %km
vsati_Galileo9 = atan2((sqrt(1-(Galileo9(3)^2)))*sin(Esati_Galileo9),(cos(Esati_Galileo9)-Galileo9(3))); %rad
rorbsati_Galileo9 = rsati_Galileo9_km .* [cos(vsati_Galileo9);sin(vsati_Galileo9);zeros(1,1440)];%km
rECI_Galileo9 = rot3d(-Galileo9(2),3)*rot3d(-Galileo9(1),1)*rot3d(-Galileo9(4),3)*rorbsati_Galileo9; %km

Mi_Beidou3(i) = t_d(i) * Beidou3(6) * (2*pi); %rad
Msati_Beidou3(i) = Beidou3(5)/180*pi + Mi_Beidou3(i); %rad
Esati_Beidou3 = kepler_e_9(Msati_Beidou3,Beidou3(3)); %rad
rsati_Beidou3_km = a_Beidou3_km*(1-Beidou3(3)*cos(Esati_Beidou3)); %km
vsati_Beidou3 = atan2((sqrt(1-(Beidou3(3)^2)))*sin(Esati_Beidou3),(cos(Esati_Beidou3)-Beidou3(3))); %rad
rorbsati_Beidou3 = rsati_Beidou3_km .* [cos(vsati_Beidou3);sin(vsati_Beidou3);zeros(1,1440)];%km
rECI_Beidou3 = rot3d(-Beidou3(2),3)*rot3d(-Beidou3(1),1)*rot3d(-Beidou3(4),3)*rorbsati_Beidou3; %km

Mi_Meteosat11(i) = t_d(i) * Meteosat11(6) * (2*pi); %rad
Msati_Meteosat11(i) = Meteosat11(5)/180*pi + Mi_Meteosat11(i); %rad
Esati_Meteosat11 = kepler_e_9(Msati_Meteosat11,Meteosat11(3)); %rad
rsati_Meteosat11_km = a_Meteosat11_km*(1-Meteosat11(3)*cos(Esati_Meteosat11)); %km
vsati_Meteosat11 = atan2((sqrt(1-(Meteosat11(3)^2)))*sin(Esati_Meteosat11),(cos(Esati_Meteosat11)-Meteosat11(3))); %rad
rorbsati_Meteosat11 = rsati_Meteosat11_km .* [cos(vsati_Meteosat11);sin(vsati_Meteosat11);zeros(1,1440)];%km
rECI_Meteosat11 = rot3d(-Meteosat11(2),3)*rot3d(-Meteosat11(1),1)*rot3d(-Meteosat11(4),3)*rorbsati_Meteosat11; %km

figure (1)
[x,y,z] = sphere(180); x=6371.*x; y=6371.*y; z=6371.*z; %km
surf(x,y,z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
hold on
plot3(rECI_ISS(1,:),rECI_ISS(2,:),rECI_ISS(3,:),'r') %km
plot3(rECI_Sentinel1B(1,:),rECI_Sentinel1B(2,:),rECI_Sentinel1B(3,:),'g') %km
plot3(rECI_Molniya1T(1,:),rECI_Molniya1T(2,:),rECI_Molniya1T(3,:),'b') %km
plot3(rECI_Galileo9(1,:),rECI_Galileo9(2,:),rECI_Galileo9(3,:),'c') %km
plot3(rECI_Beidou3(1,:),rECI_Beidou3(2,:),rECI_Beidou3(3,:),'m') %km
plot3(rECI_Meteosat11(1,:),rECI_Meteosat11(2,:),rECI_Meteosat11(3,:),'k') %km
title('24 h ECI orbit')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('','ISS','Sentinel 1B','Molniya-1T','Galileo 9','Beidou IGSO 3','Meteosat 11')
hold off

%ECI velocity

p_ofellipse_ISS = rsati_ISS_km .* (1 + ISS(3)*cos(vsati_ISS)); %km
Mag_C_ISS = sqrt(p_ofellipse_ISS*GM);%km^2/s
Unitvec_C_ISS = [sind(ISS(2))*sind(ISS(1)),-cosd(ISS(2))*sind(ISS(1)),cosd(ISS(1))]';
Unitvec_K_ISS = [cosd(ISS(2)),sind(ISS(2)),0]';
P_ISS = cosd(ISS(4))*Unitvec_K_ISS + sind(ISS(4))*cross(Unitvec_C_ISS,Unitvec_K_ISS);
Q_ISS = -sind(ISS(4))*Unitvec_K_ISS + cosd(ISS(4))*cross(Unitvec_C_ISS,Unitvec_K_ISS);
vorb_ISS = (Mag_C_ISS/p_ofellipse_ISS)*(-sin(vsati_ISS).*P_ISS + ((ISS(3)+cos(vsati_ISS)).*Q_ISS));
for i = 1:1440
    v_mag_ISS(1,i) = sqrt(vorb_ISS(1,i)^2+vorb_ISS(2,i)^2+vorb_ISS(3,i)^2);
end

p_ofellipse_Sentinel1B = rsati_Sentinel1B_km .* (1 + Sentinel1B(3)*cos(vsati_Sentinel1B)); %km
Mag_C_Sentinel1B = sqrt(p_ofellipse_Sentinel1B*GM);%km^2/s
Unitvec_C_Sentinel1B = [sind(Sentinel1B(2))*sind(Sentinel1B(1)),-cosd(Sentinel1B(2))*sind(Sentinel1B(1)),cosd(Sentinel1B(1))]';
Unitvec_K_Sentinel1B = [cosd(Sentinel1B(2)),sind(Sentinel1B(2)),0]';
P_Sentinel1B = cosd(Sentinel1B(4))*Unitvec_K_Sentinel1B + sind(Sentinel1B(4))*cross(Unitvec_C_Sentinel1B,Unitvec_K_Sentinel1B);
Q_Sentinel1B = -sind(Sentinel1B(4))*Unitvec_K_Sentinel1B + cosd(Sentinel1B(4))*cross(Unitvec_C_Sentinel1B,Unitvec_K_Sentinel1B);
vorb_Sentinel1B = (Mag_C_Sentinel1B/p_ofellipse_Sentinel1B)*(-sin(vsati_Sentinel1B).*P_Sentinel1B + ((Sentinel1B(3)+cos(vsati_Sentinel1B)).*Q_Sentinel1B));
for i = 1:1440
    v_mag_Sentinel1B(1,i) = sqrt(vorb_Sentinel1B(1,i)^2+vorb_Sentinel1B(2,i)^2+vorb_Sentinel1B(3,i)^2);
end

p_ofellipse_Molniya1T = rsati_Molniya1T_km .* (1 + Molniya1T(3)*cos(vsati_Molniya1T)); %km
Mag_C_Molniya1T = sqrt(p_ofellipse_Molniya1T*GM);%km^2/s
Unitvec_C_Molniya1T = [sind(Molniya1T(2))*sind(Molniya1T(1)),-cosd(Molniya1T(2))*sind(Molniya1T(1)),cosd(Molniya1T(1))]';
Unitvec_K_Molniya1T = [cosd(Molniya1T(2)),sind(Molniya1T(2)),0]';
P_Molniya1T = cosd(Molniya1T(4))*Unitvec_K_Molniya1T + sind(Molniya1T(4))*cross(Unitvec_C_Molniya1T,Unitvec_K_Molniya1T);
Q_Molniya1T = -sind(Molniya1T(4))*Unitvec_K_Molniya1T + cosd(Molniya1T(4))*cross(Unitvec_C_Molniya1T,Unitvec_K_Molniya1T);
vorb_Molniya1T = (Mag_C_Molniya1T/p_ofellipse_Molniya1T)*(-sin(vsati_Molniya1T).*P_Molniya1T + ((Molniya1T(3)+cos(vsati_Molniya1T)).*Q_Molniya1T));
for i = 1:1440
    v_mag_Molniya1T(1,i) = sqrt(vorb_Molniya1T(1,i)^2+vorb_Molniya1T(2,i)^2+vorb_Molniya1T(3,i)^2);
end

p_ofellipse_Galileo9 = rsati_Galileo9_km .* (1 + Galileo9(3)*cos(vsati_Galileo9)); %km
Mag_C_Galileo9 = sqrt(p_ofellipse_Galileo9*GM);%km^2/s
Unitvec_C_Galileo9 = [sind(Galileo9(2))*sind(Galileo9(1)),-cosd(Galileo9(2))*sind(Galileo9(1)),cosd(Galileo9(1))]';
Unitvec_K_Galileo9 = [cosd(Galileo9(2)),sind(Galileo9(2)),0]';
P_Galileo9 = cosd(Galileo9(4))*Unitvec_K_Galileo9 + sind(Galileo9(4))*cross(Unitvec_C_Galileo9,Unitvec_K_Galileo9);
Q_Galileo9 = -sind(Galileo9(4))*Unitvec_K_Galileo9 + cosd(Galileo9(4))*cross(Unitvec_C_Galileo9,Unitvec_K_Galileo9);
vorb_Galileo9 = (Mag_C_Galileo9/p_ofellipse_Galileo9)*(-sin(vsati_Galileo9).*P_Galileo9 + ((Galileo9(3)+cos(vsati_Galileo9)).*Q_Galileo9));
for i = 1:1440
    v_mag_Galileo9(1,i) = sqrt(vorb_Galileo9(1,i)^2+vorb_Galileo9(2,i)^2+vorb_Galileo9(3,i)^2);
end

p_ofellipse_Beidou3 = rsati_Beidou3_km .* (1 + Beidou3(3)*cos(vsati_Beidou3)); %km
Mag_C_Beidou3 = sqrt(p_ofellipse_Beidou3*GM);%km^2/s
Unitvec_C_Beidou3 = [sind(Beidou3(2))*sind(Beidou3(1)),-cosd(Beidou3(2))*sind(Beidou3(1)),cosd(Beidou3(1))]';
Unitvec_K_Beidou3 = [cosd(Beidou3(2)),sind(Beidou3(2)),0]';
P_Beidou3 = cosd(Beidou3(4))*Unitvec_K_Beidou3 + sind(Beidou3(4))*cross(Unitvec_C_Beidou3,Unitvec_K_Beidou3);
Q_Beidou3 = -sind(Beidou3(4))*Unitvec_K_Beidou3 + cosd(Beidou3(4))*cross(Unitvec_C_Beidou3,Unitvec_K_Beidou3);
vorb_Beidou3 = (Mag_C_Beidou3/p_ofellipse_Beidou3)*(-sin(vsati_Beidou3).*P_Beidou3 + ((Beidou3(3)+cos(vsati_Beidou3)).*Q_Beidou3));
for i = 1:1440
    v_mag_Beidou3(1,i) = sqrt(vorb_Beidou3(1,i)^2+vorb_Beidou3(2,i)^2+vorb_Beidou3(3,i)^2);
end

p_ofellipse_Meteosat11 = rsati_Meteosat11_km .* (1 + Meteosat11(3)*cos(vsati_Meteosat11)); %km
Mag_C_Meteosat11 = sqrt(p_ofellipse_Meteosat11*GM);%km^2/s
Unitvec_C_Meteosat11 = [sind(Meteosat11(2))*sind(Meteosat11(1)),-cosd(Meteosat11(2))*sind(Meteosat11(1)),cosd(Meteosat11(1))]';
Unitvec_K_Meteosat11 = [cosd(Meteosat11(2)),sind(Meteosat11(2)),0]';
P_Meteosat11 = cosd(Meteosat11(4))*Unitvec_K_Meteosat11 + sind(Meteosat11(4))*cross(Unitvec_C_Meteosat11,Unitvec_K_Meteosat11);
Q_Meteosat11 = -sind(Meteosat11(4))*Unitvec_K_Meteosat11 + cosd(Meteosat11(4))*cross(Unitvec_C_Meteosat11,Unitvec_K_Meteosat11);
vorb_Meteosat11 = (Mag_C_Meteosat11/p_ofellipse_Meteosat11)*(-sin(vsati_Meteosat11).*P_Meteosat11 + ((Meteosat11(3)+cos(vsati_Meteosat11)).*Q_Meteosat11));
for i = 1:1440
    v_mag_Meteosat11(1,i) = sqrt(vorb_Meteosat11(1,i)^2+vorb_Meteosat11(2,i)^2+vorb_Meteosat11(3,i)^2);
end

figure (2)
xaxis = linspace(0,24,1440);
plot(xaxis,vorb_ISS(1,:))
hold on 
plot(xaxis,vorb_ISS(2,:))
plot(xaxis,vorb_ISS(3,:))
title('24 h ECI velocity components of ISS')
xlabel('time [h]')
ylabel('velocity [km/s]')
legend('vel x','vel y','vel z')
hold off
figure (3)
plot(xaxis,v_mag_ISS)
xlabel('time [h]')
ylabel('velocity [km/s]')
title('24 h ECI velocity magintude of ISS')


figure (4)
xaxis = linspace(0,24,1440);
plot(xaxis,vorb_Sentinel1B(1,:),'r')
hold on 
plot(xaxis,vorb_Sentinel1B(2,:),'g')
plot(xaxis,vorb_Sentinel1B(3,:),'b')
title('24 h ECI velocity components of Sentinel 1B')
xlabel('time [h]')
ylabel('velocity [km/s]')
legend('vel x','vel y','vel z')
hold off
figure (5)
plot(xaxis,v_mag_Sentinel1B)
xlabel('time [h]')
ylabel('velocity [km/s]')
title('24 h ECI velocity magintude of Sentinel 1B')

figure (6)
xaxis = linspace(0,24,1440);
plot(xaxis,vorb_Molniya1T(1,:))
hold on 
plot(xaxis,vorb_Molniya1T(2,:))
plot(xaxis,vorb_Molniya1T(3,:))
title('24 h ECI velocity components of Molniya-1T')
xlabel('time [h]')
ylabel('velocity [km/s]')
legend('vel x','vel y','vel z')
hold off
figure (7)
plot(xaxis,v_mag_Molniya1T)
xlabel('time [h]')
ylabel('velocity [km/s]')
title('24 h ECI velocity magintude of Molniya-1T')

figure (8)
xaxis = linspace(0,24,1440);
plot(xaxis,vorb_Galileo9(1,:))
hold on 
plot(xaxis,vorb_Galileo9(2,:))
plot(xaxis,vorb_Galileo9(3,:))
title('24 h ECI velocity components of Galileo 9')
xlabel('time [h]')
ylabel('velocity [km/s]')
legend('vel x','vel y','vel z')
hold off
figure (9)
plot(xaxis,v_mag_Galileo9)
xlabel('time [h]')
ylabel('velocity [km/s]')
title('24 h ECI velocity magintude of Galileo 9')


figure (10)
xaxis = linspace(0,24,1440);
plot(xaxis,vorb_Beidou3(1,:))
hold on 
plot(xaxis,vorb_Beidou3(2,:))
plot(xaxis,vorb_Beidou3(3,:))
title('24 h ECI velocity components of Beidou IGSO 3')
xlabel('time [h]')
ylabel('velocity [km/s]')
legend('vel x','vel y','vel z')
hold off
figure (11)
plot(xaxis,v_mag_Beidou3)
xlabel('time [h]')
ylabel('velocity [km/s]')
title('24 h ECI velocity magintude of Beidou IGSO 3')

figure (12)
xaxis = linspace(0,24,1440);
plot(xaxis,vorb_Meteosat11(1,:))
hold on 
plot(xaxis,vorb_Meteosat11(2,:))
plot(xaxis,vorb_Meteosat11(3,:))
title('24 h ECI velocity components of Meteosat 11')
xlabel('time [h]')
ylabel('velocity [km/s]')
legend('vel x','vel y','vel z')
hold off
figure (13)
plot(xaxis,v_mag_Meteosat11)
xlabel('time [h]')
ylabel('velocity [km/s]')
title('24 h ECI velocity magintude of Meteosat 11')


%% 4

%incresing 1 min to reach 24 hours
j = 1:1440;
t_ISS(j) = t_ISS + minutes(j-1); 
t_Sentinel1B(j) = t_Sentinel1B + minutes(j-1); 
t_Molniya1T(j) = t_Molniya1T + minutes(j-1); 
t_Galileo9(j) = t_Galileo9 + minutes(j-1); 
t_Beidou3(j) = t_Beidou3 + minutes(j-1); 
t_Meteosat11(j) = t_Meteosat11 + minutes(j-1); 

%[JD_UT1,MJD_UT1] = gre2jd(year(UT1),month(UT1),day(UT1),hour(UT1),minute(UT1),second(UT1));
for j = 1:1440
    JD_UT1_ISS = gre2jd(year(t_ISS),month(t_ISS),day(t_ISS),hour(t_ISS),minute(t_ISS),second(t_ISS));
    JD_UT1_Sentinel1B = gre2jd(year(t_Sentinel1B),month(t_Sentinel1B),day(t_Sentinel1B),hour(t_Sentinel1B),minute(t_Sentinel1B),second(t_Sentinel1B));
    JD_UT1_Molniya1T = gre2jd(year(t_Molniya1T),month(t_Molniya1T),day(t_Molniya1T),hour(t_Molniya1T),minute(t_Molniya1T),second(t_Molniya1T));
    JD_UT1_Galileo9 = gre2jd(year(t_Galileo9),month(t_Galileo9),day(t_Galileo9),hour(t_Galileo9),minute(t_Galileo9),second(t_Galileo9));
    JD_UT1_Beidou3 = gre2jd(year(t_Beidou3),month(t_Beidou3),day(t_Beidou3),hour(t_Beidou3),minute(t_Beidou3),second(t_Beidou3));
    JD_UT1_Meteosat11 = gre2jd(year(t_Meteosat11),month(t_Meteosat11),day(t_Meteosat11),hour(t_Meteosat11),minute(t_Meteosat11),second(t_Meteosat11));
end
JD_UT1_ISS = JD_UT1_ISS(1:1440);
JD_UT1_Sentinel1B = JD_UT1_Sentinel1B(1:1440);
JD_UT1_Molniya1T = JD_UT1_Molniya1T(1:1440);
JD_UT1_Galileo9 = JD_UT1_Galileo9(1:1440);
JD_UT1_Beidou3 = JD_UT1_Beidou3(1:1440);
JD_UT1_Meteosat11 = JD_UT1_Meteosat11(1:1440);

%GMST 

tt_ISS = (JD_UT1_ISS - 2451545.0)/36525; 
GMST_ISS = (F(JD_UT1_ISS)*86400 + 24110.54841 - 86400/2 + ...
        8640184.812866 .* tt_ISS + 0.093104 .* tt_ISS .* tt_ISS - ...
        6.2e-6 .* tt_ISS .* tt_ISS .* tt_ISS)/3600; % hr
GMST_ISS = rem(GMST_ISS,24);% hr decimal degree
hdGMST_ISS = GMST_ISS*15; % hr angle degree

tt_Sentinel1B = (JD_UT1_Sentinel1B - 2451545.0)/36525; 
GMST_Sentinel1B = (F(JD_UT1_Sentinel1B)*86400 + 24110.54841 - 86400/2 + ...
        8640184.812866 .* tt_Sentinel1B + 0.093104 .* tt_Sentinel1B .* tt_Sentinel1B - ...
        6.2e-6 .* tt_Sentinel1B .* tt_Sentinel1B .* tt_Sentinel1B)/3600; % hr
GMST_Sentinel1B = rem(GMST_Sentinel1B,24);% hr decimal degree
hdGMST_Sentinel1B = GMST_Sentinel1B*15; % hr angle degree

tt_Molniya1T = (JD_UT1_Molniya1T - 2451545.0)/36525; 
GMST_Molniya1T = (F(JD_UT1_Molniya1T)*86400 + 24110.54841 - 86400/2 + ...
        8640184.812866 .* tt_Molniya1T + 0.093104 .* tt_Molniya1T .* tt_Molniya1T - ...
        6.2e-6 .* tt_Molniya1T .* tt_Molniya1T .* tt_Molniya1T)/3600; % hr
GMST_Molniya1T = rem(GMST_Molniya1T,24);% hr decimal degree
hdGMST_Molniya1T = GMST_Molniya1T*15; % hr angle degree

tt_Galileo9 = (JD_UT1_Galileo9 - 2451545.0)/36525; 
GMST_Galileo9 = (F(JD_UT1_Galileo9)*86400 + 24110.54841 - 86400/2 + ...
        8640184.812866 .* tt_Galileo9 + 0.093104 .* tt_Galileo9 .* tt_Galileo9 - ...
        6.2e-6 .* tt_Galileo9 .* tt_Galileo9 .* tt_Galileo9)/3600; % hr
GMST_Galileo9 = rem(GMST_Galileo9,24);% hr decimal degree
hdGMST_Galileo9 = GMST_Galileo9*15; % hr angle degree

tt_Beidou3 = (JD_UT1_Beidou3 - 2451545.0)/36525; 
GMST_Beidou3 = (F(JD_UT1_Beidou3)*86400 + 24110.54841 - 86400/2 + ...
        8640184.812866 .* tt_Beidou3 + 0.093104 .* tt_Beidou3 .* tt_Beidou3 - ...
        6.2e-6 .* tt_Beidou3 .* tt_Beidou3 .* tt_Beidou3)/3600; % hr
GMST_Beidou3 = rem(GMST_Beidou3,24);% hr decimal degree
hdGMST_Beidou3 = GMST_Beidou3*15; % hr angle degree

tt_Meteosat11 = (JD_UT1_Meteosat11 - 2451545.0)/36525; 
GMST_Meteosat11 = (F(JD_UT1_Meteosat11)*86400 + 24110.54841 - 86400/2 + ...
        8640184.812866 .* tt_Meteosat11 + 0.093104 .* tt_Meteosat11 .* tt_Meteosat11 - ...
        6.2e-6 .* tt_Meteosat11 .* tt_Meteosat11 .* tt_Meteosat11)/3600; % hr
GMST_Meteosat11 = rem(GMST_Meteosat11,24);% hr decimal degree
hdGMST_Meteosat11 = GMST_Meteosat11*15; % hr angle degree

%Earth Rototation
for i = 1:1440
    
    R_ISS = rot3d(hdGMST_ISS(i),3);
    rECEF_ISS(1:3,i) = R_ISS*(rECI_ISS(1:3,i)); %km
    R_Sentinel1B  = rot3d(hdGMST_Sentinel1B(i),3);
    rECEF_Sentinel1B(1:3,i) = R_Sentinel1B*(rECI_Sentinel1B(1:3,i)); %km
    R_Molniya1T = rot3d(hdGMST_Molniya1T(i),3);
    rECEF_Molniya1T(1:3,i) = R_Molniya1T*(rECI_Molniya1T(1:3,i)); %km
    R_Galileo9 = rot3d(hdGMST_Galileo9(i),3);
    rECEF_Galileo9(1:3,i) = R_Galileo9*(rECI_Galileo9(1:3,i)); %km
    R_Beidou3 = rot3d(hdGMST_Beidou3(i),3);
    rECEF_Beidou3(1:3,i) = R_Beidou3*(rECI_Beidou3(1:3,i)); %km
    R_Meteosat11 = rot3d(hdGMST_Meteosat11(i),3);
    rECEF_Meteosat11(1:3,i) = R_Meteosat11*(rECI_Meteosat11(1:3,i)); %km
end

figure (14)
[x,y,z] = sphere(180); x=6371.*x; y=6371.*y; z=6371.*z; %km
earth=surf(x,y,-z,'FaceColor','none','EdgeColor',[1 1 1]);
cdata = imread('my_Earth_image.jpg');
set(earth,'FaceColor','texturemap','CData',cdata,...
'FaceAlpha',1,'EdgeColor','none');
axis equal
axis auto
axis vis3d
hold on
plot3(rECEF_ISS(1,:),rECEF_ISS(2,:),rECEF_ISS(3,:),'r') %km
plot3(rECEF_Sentinel1B(1,:),rECEF_Sentinel1B(2,:),rECEF_Sentinel1B(3,:),'g') %km
plot3(rECEF_Molniya1T(1,:),rECEF_Molniya1T(2,:),rECEF_Molniya1T(3,:),'b') %km
plot3(rECEF_Galileo9(1,:),rECEF_Galileo9(2,:),rECEF_Galileo9(3,:),'c') %km
plot3(rECEF_Beidou3(1,:),rECEF_Beidou3(2,:),rECEF_Beidou3(3,:),'m') %km
plot3(rECEF_Meteosat11(1,:),rECEF_Meteosat11(2,:),rECEF_Meteosat11(3,:),'k') %km
title('24 h ECEF orbit')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('','ISS','Sentinel 1B','Molniya-1T','Galileo 9','Beidou IGSO 3','Meteosat 11')
hold off

%% 5 
for i = 1:1440
    
    plr_ISS_1 = xyz2latlonr([rECEF_ISS(1,i),rECEF_ISS(2,i),rECEF_ISS(3,i)]);
    plr_ISS_1 = plr_ISS_1(1:2); % get phi, lam [rad]
    pl_ISS_all(1:2,i) = plr_ISS_1'; % [rad]
    
    plr_Sentinel1B_1 = xyz2latlonr([rECEF_Sentinel1B(1,i),rECEF_Sentinel1B(2,i),rECEF_Sentinel1B(3,i)]);
    plr_Sentinel1B_1 = plr_Sentinel1B_1(1:2); % get phi, lam [rad]
    pl_Sentinel1B_all(1:2,i) = plr_Sentinel1B_1'; % [rad]

    plr_Molniya1T_1 = xyz2latlonr([rECEF_Molniya1T(1,i),rECEF_Molniya1T(2,i),rECEF_Molniya1T(3,i)]);
    plr_Molniya1T_1 = plr_Molniya1T_1(1:2); % get phi, lam [rad]
    pl_Molniya1T_all(1:2,i) = plr_Molniya1T_1'; % [rad]

    plr_Galileo9_1 = xyz2latlonr([rECEF_Galileo9(1,i),rECEF_Galileo9(2,i),rECEF_Galileo9(3,i)]);
    plr_Galileo9_1 = plr_Galileo9_1(1:2); % get phi, lam [rad]
    pl_Galileo9_all(1:2,i) = plr_Galileo9_1'; % [rad]

    plr_Beidou3_1 = xyz2latlonr([rECEF_Beidou3(1,i),rECEF_Beidou3(2,i),rECEF_Beidou3(3,i)]);
    plr_Beidou3_1 = plr_Beidou3_1(1:2); % get phi, lam [rad]
    pl_Beidou3_all(1:2,i) = plr_Beidou3_1'; % [rad]

    plr_Meteosat11_1 = xyz2latlonr([rECEF_Meteosat11(1,i),rECEF_Meteosat11(2,i),rECEF_Meteosat11(3,i)]);
    plr_Meteosat11_1 = plr_Meteosat11_1(1:2); % get phi, lam [rad]
    pl_Meteosat11_all(1:2,i) = plr_Meteosat11_1'; % [rad]    
end

% rad -> deg
% phi,lam
pl_ISS_all_deg = pl_ISS_all*(180/pi); 
pl_Sentinel1B_all_deg = pl_Sentinel1B_all*(180/pi); 
pl_Molniya1T_all_deg = pl_Molniya1T_all*(180/pi);
pl_Galileo9_all_deg = pl_Galileo9_all*(180/pi);
pl_Beidou3_all_deg = pl_Beidou3_all*(180/pi);
pl_Meteosat11_all_deg = pl_Meteosat11_all*(180/pi);

figure (15)
cdata = imread('my_Earth_image.jpg');
imagesc([-180 180], [90 -90], cdata);
set(gca,'ydir','normal');
hold on
plot(pl_ISS_all_deg(2,:),pl_ISS_all_deg(1,:),'.r')
plot(pl_Sentinel1B_all_deg(2,:),pl_Sentinel1B_all_deg(1,:),'.g')
plot(pl_Molniya1T_all_deg(2,:),pl_Molniya1T_all_deg(1,:),'.b')
plot(pl_Galileo9_all_deg(2,:),pl_Galileo9_all_deg(1,:),'.c')
plot(pl_Beidou3_all_deg(2,:),pl_Beidou3_all_deg(1,:),'.m')
plot(pl_Meteosat11_all_deg(2,:),pl_Meteosat11_all_deg(1,:),'.y')
title('24 ground track six Earth satellites')
xlabel('longitude [deg]')
ylabel('latitude [deg]')
legend('ISS','Sentinel 1B','Molniya-1T','Galileo 9','Beidou IGSO 3','Meteosat 11','Location','northeastoutside')
hold off

%% 6

% Cartesian Berlin [m.] -> Spherical Berlin [deg]
Berlin_cat_m =[3782971.01,902152.49,5038375.59]; % [m.]
Berlin_plr_rad = xyz2latlonr(Berlin_cat_m); %[rad.]
Berlin_pl_deg = Berlin_plr_rad(1:2) * (180/pi); %[deg.]

% ECEF -> Local Horizontal System
for i = 1:1440
    
    rtopo_ISS(1:3,i) = ref3d(1)*rot3d(90-Berlin_pl_deg(1),2)*rot3d(Berlin_pl_deg(2),3)*rECEF_ISS(1:3,i); %[km]
    z_ISS(i) = (pi/2)-atan(rtopo_ISS(3,i)/sqrt(rtopo_ISS(1,i)^2+rtopo_ISS(2,i)^2));%rad
    dist_3D_ISS(i) = sqrt(rtopo_ISS(1,i)^2+rtopo_ISS(2,i)^2+rtopo_ISS(3,i)^2) * 1000; %m
    
    rtopo_Sentinel1B(1:3,i) = ref3d(1)*rot3d(90-Berlin_pl_deg(1),2)*rot3d(Berlin_pl_deg(2),3)*rECEF_Sentinel1B(1:3,i); %[km]
    z_Sentinel1B(i) = (pi/2)-atan(rtopo_Sentinel1B(3,i)/sqrt(rtopo_Sentinel1B(1,i)^2+rtopo_Sentinel1B(2,i)^2));%rad
    dist_3D_Sentinel1B(i) = sqrt(rtopo_Sentinel1B(1,i)^2+rtopo_Sentinel1B(2,i)^2+rtopo_Sentinel1B(3,i)^2) * 1000; %m

    rtopo_Molniya1T(1:3,i) = ref3d(1)*rot3d(90-Berlin_pl_deg(1),2)*rot3d(Berlin_pl_deg(2),3)*rECEF_Molniya1T(1:3,i); %[km]
    z_Molniya1T(i) = (pi/2)-atan(rtopo_Molniya1T(3,i)/sqrt(rtopo_Molniya1T(1,i)^2+rtopo_Molniya1T(2,i)^2));%rad
    dist_3D_Molniya1T(i) = sqrt(rtopo_Molniya1T(1,i)^2+rtopo_Molniya1T(2,i)^2+rtopo_Molniya1T(3,i)^2) * 1000; %m

    rtopo_Galileo9(1:3,i) = ref3d(1)*rot3d(90-Berlin_pl_deg(1),2)*rot3d(Berlin_pl_deg(2),3)*rECEF_Galileo9(1:3,i); %[km]
    z_Galileo9(i) = (pi/2)-atan(rtopo_Galileo9(3,i)/sqrt(rtopo_Galileo9(1,i)^2+rtopo_Galileo9(2,i)^2));%rad
    dist_3D_Galileo9(i) = sqrt(rtopo_Galileo9(1,i)^2+rtopo_Galileo9(2,i)^2+rtopo_Galileo9(3,i)^2) * 1000; %m

    rtopo_Beidou3(1:3,i) = ref3d(1)*rot3d(90-Berlin_pl_deg(1),2)*rot3d(Berlin_pl_deg(2),3)*rECEF_Beidou3(1:3,i); %[km]
    z_Beidou3(i) = (pi/2)-atan(rtopo_Beidou3(3,i)/sqrt(rtopo_Beidou3(1,i)^2+rtopo_Beidou3(2,i)^2));%rad
    dist_3D_Beidou3(i) = sqrt(rtopo_Beidou3(1,i)^2+rtopo_Beidou3(2,i)^2+rtopo_Beidou3(3,i)^2) * 1000; %m

    rtopo_Meteosat11(1:3,i) = ref3d(1)*rot3d(90-Berlin_pl_deg(1),2)*rot3d(Berlin_pl_deg(2),3)*rECEF_Meteosat11(1:3,i); %[km]
    z_Meteosat11(i) = (pi/2)-atan(rtopo_Meteosat11(3,i)/sqrt(rtopo_Meteosat11(1,i)^2+rtopo_Meteosat11(2,i)^2));%rad
    dist_3D_Meteosat11(i) = sqrt(rtopo_Meteosat11(1,i)^2+rtopo_Meteosat11(2,i)^2+rtopo_Meteosat11(3,i)^2) * 1000; %m

end

min_z_ISS = min(z_ISS) %rad
epoch_min_z_ISS = find(z_ISS <= min_z_ISS) %index

min_dist_3D_ISS = min(dist_3D_ISS) %m
epoch_min_dist3D_ISS = find(dist_3D_ISS <= min_dist_3D_ISS) %index

