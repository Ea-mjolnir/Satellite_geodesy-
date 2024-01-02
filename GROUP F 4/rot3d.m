function [R_i] = rot3d(angle,axis)
%% ROT3D creates a [3x3] rotation matrix
%    usage:
%    R_i = rot3d(angle,axis)
%    input:
%    angle rotation angle [deg.dec]
%    axis  [3x3] rotation matrix
%    output:
%    R_i   [3x3] rotation matrix
%    external calls:
%    none
%% robustness
if axis~=1 && axis~=2 && axis~=3
   error('2nd input argument "axis" must be 1, 2, or 3.');
end
if ~isscalar(angle)
   error('1st input argument "angle" must be scalar.');
end


%% computations
R_i = eye(3);
sang = sind(angle);
cang = cosd(angle);
switch axis
    case 1
       R_i(2,2) = cang;
       R_i(2,3) = sang;
       R_i(3,2) = -sang;
       R_i(3,3) = cang;
    case 2
       R_i(1,1) = cang;
       R_i(1,3) = -sang;
       R_i(3,1) = sang;
       R_i(3,3) = cang;
    case 3
       R_i(1,1) = cang;
       R_i(1,2) = sang;
       R_i(2,1) = -sang;
       R_i(2,2) = cang;
end
