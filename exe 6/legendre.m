function Pnm = legendre(n,m,theta)

% input degree n, order m, and colatitude theta (deg)
% checks and stuff

if ~isscalar(n) || ~isscalar(m)
    error('n and m have to be integer scalars! aborting');
end
if ~isequal(size(theta,1),1) && ~isequal(size(theta,2),1)
    error('\vartheta has to be a vector, not a matrix!');
end
if m>n, error('m > n cant be!'); 
end

a1 = sqrt(3);
an = sqrt((2*n+1)/(2*n));
bnm = sqrt(((2*n+1)*(2*n-1))/((n+m)*(n-m)));
cnm = sqrt(((2*n+1)*(n-m-1)*(n+m-1))/((2*n-3)*(n+m)*(n-m)));

if n == 0 & m == 0
    Pnm = 1;
elseif n == n & m == n
    Pnm = an*sind(theta)*legendre(n-1,n-1,theta);
elseif n == n & m == (n-1)
    Pnm = bnm*cosd(theta)*legendre(n-1,n-1,theta);
else
    Pnm = bnm*cosd(theta)*legendre(n-1,m,theta) - cnm*legendre(n-2,m,theta);
end