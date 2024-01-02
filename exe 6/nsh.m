function CnmSnm = nsh(n,m,theta,lambda)
% usage [Cnm,Snm] = NSH(n,m,theta,lambda)
% NSH computes normalized spherical harmonics with
% input degree n, order m, colatitude theta (deg), longitude lambda (deg)
% input theta and lambda can be in vectorial form
% n and m must be scalar, with n>=m
if m > n, error('Check input! m must be smaller or equal n!'); 
end
Pnm = legendre(n,m,theta);
Cnm = Pnm*cosd(m*lambda);
Snm = Pnm*sind(m*lambda);
CnmSnm = [Cnm,Snm];