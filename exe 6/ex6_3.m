%% Ex 3
clc
clear

n = [2,2,2,3,3,3,3];
m = [0,1,2,0,1,2,3];
lambda=0:2:360; %deg
theta=0:1:180;  %deg

GM = 398600.44; %[km^3/s^2]
R = 6378; %[Km]
r = 6378; %[km]
dCnm = 10^-6 * [0.002,0,2.439,0.957,2.029,0.904,0.721];
dSnm = 10^-6 * [10^6,0,-1.400,10^6,0.249,-0.619,1.144];
for i = 1:length(theta)
    for j = 1:7
     Pnm(j,i) = legendre(n(j),m(j),theta(i));
     T_dispo(j,i) = (GM/r) * sum((R/r)^3+(R/r)^3+(R/r)^3+(R/r)^4+(R/r)^4+ ...
        (R/r)^4+(R/r)^4) * ...
        sum(Pnm(j,i)*(dCnm(j)*cosd(m(j)*lambda(i))+dSnm(j)*sind(m(j)*lambda(i))));
    end
end
