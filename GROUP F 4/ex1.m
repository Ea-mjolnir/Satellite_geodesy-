clc
clear
% vernaleqinox
alpha = 0; delta = 0;
vernaleqinox = [delta,alpha,1];
xi0 = latlonr2xyz(vernaleqinox);  %(rad,rad,m)-->(m,m,m)

% JD in 1990-2030
yyyy = 1990:2030;
mm=1; dd=1; ut1=12; minute=0; second=0;
jd = gre2jd(yyyy,mm,dd,ut1,minute,second);
jd = jd(1:41); % Extract only Julian Date

% P(JD)
%xip = 
for i = 1:length(jd)
    Pjd = prec(jd(i)); 
    xip = Pjd * xi0';
    % (m,m,m) --> (rad,rad,m)
    sph_ip = xyz2latlonr(xip); 
    % [delta,alpha,1] - [rad,rad,m]
    
    delta2(i) = sph_ip(1)*180/pi; % degree
    alpha2(i) = sph_ip(2)*12/pi; % hour angle
end
delta2;
alpha2;
yyyy;

fig1 = plotyy(yyyy,alpha2,yyyy,delta2);
title('precession of a fictitious celestial object at (\alpha=0;\delta=0)')
xlabel('time (year)')
ylabel(fig1(1),'\alpha (h)')
ylabel(fig1(2),'\delta (deg)')