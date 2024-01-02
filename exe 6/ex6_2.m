%% Ex 2
clc
clear
n = [4,5,5];
m = [0,5,4];
lambda=0:2:360; %deg
theta=0:1:180;  %deg
for i = 1:length(theta)
    for j = 1:3
    CnmSnm = nsh(n(j),m(j),theta(i),lambda(i));
    Cnm(j,i) = CnmSnm(1);
    Snm(j,i) = CnmSnm(2);
    end
end

%% Zonal n=4 m=0
figure(1)
imagesc(120:240,theta,Cnm(1,:)')
colorbar
title('Zonal: Cnm, n = 4, m = 0')
xlabel('longitude [deg]')
ylabel('co-latitude [deg]')

figure(2)
imagesc(120:240,theta,Snm(1,:)')
colorbar
title('Zonal: Snm, n = 4, m = 0')
xlabel('longitude [deg]')
ylabel('co-latitude [deg]')

%% Sectorial n=5 m=5
figure(3)
imagesc(120:240,theta,Cnm(2,:)')
colorbar
title('Sectorial: Cnm, n = 5, m = 5')
xlabel('longitude [deg]')
ylabel('co-latitude [deg]')

figure(4)
imagesc(120:240,theta,Snm(2,:)')
colorbar
title('Sectorial: Snm, n = 5, m = 5')
xlabel('longitude [deg]')
ylabel('co-latitude [deg]')
%% Tesseral n=5 m=4
figure(5)
imagesc(120:240,theta,Cnm(3,:)')
colorbar
title('Tesseral: Cnm, n = 5, m = 4')
xlabel('longitude [deg]')
ylabel('co-latitude [deg]')

figure(6)
imagesc(120:240,theta,Snm(3,:)')
colorbar
title('Tesseral: Snm, n = 5, m = 4')
xlabel('longitude [deg]')
ylabel('co-latitude [deg]')

% get a unit sphere in correct dimensions
[x,y,z]=sphere(length(Cnm(1,:)));
% inflate sphere to Earth radius [km]
R = 6371; %[km]
x=R*x; y=R*y; z=R*z; 
% flip z-direction (MatLab specific issue)
z=-z;

figure
%surf(x,y,z,Snm')
surf(x,y,z,Cnm(1,:))
colorbar