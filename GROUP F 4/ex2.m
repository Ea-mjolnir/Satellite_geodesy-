clc
clear
% vernaleqinox
alpha = 0; delta = 0;
vernaleqinox = [delta,alpha,1];
xi0 = latlonr2xyz(vernaleqinox);  %(rad,rad,m)-->(m,m,m)

% JD in 1990-2030 with monthly resolution
mm=1:12;
dd=1; ut1=12; minute=0; second=0;
i = 1;
for yyyy = 1990:2030
    % Seperated Each Specific Year
    jd = gre2jd(yyyy,mm,dd,ut1,minute,second); 
    % JD of 12 month in Specific year
    jd = jd(1:12); % Extract only Julian Date
    jd_all(i:i+11) = jd; %Insert Each Year of JD value in 12 months
    i = i +12;
end
jd_all; % All JD in 1990-2030 with monthly resolution

%N(JD)
%xin =
for j = 1:length(jd_all)
    Njd = nut(jd_all(j)); 
    xin = Njd * xi0';
    % (m,m,m) --> (rad,rad,m)
    sph_in = xyz2latlonr(xin); 
    % [delta,alpha,1] - [rad,rad,m]
    
    delta2(j) = sph_in(1)*180/pi; % degree
    alpha2(j) = sph_in(2)*12/pi; % hour angle
end
delta2;
alpha2;
xaxis = linspace(1990,2030,492);
fig1 = plotyy(xaxis,alpha2,xaxis,delta2);
title('nutation of a fictitious celestial object at (\alpha=0;\delta=0)')
xlabel('time (year)')
ylabel(fig1(1),'\alpha (h)')
ylabel(fig1(2),'\delta (deg)')