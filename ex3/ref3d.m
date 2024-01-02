function [R_j] = ref3d(axis)
%% REF3D creates a [3x3] reflection matrix
%    usage:
%    R_j = ref3d(axis)
%    input:
%    axis  [3x3] reflection matrix
%    output:
%    R_j   [3x3] reflection matrix
%    external calls:
%    none
%% robustness
if axis~=1 && axis~=2 && axis~=3
   error('1st input argument "axis" must be 1, 2, or 3.');
end


%% computations
R_j = zeros(3);
switch axis
    case 1
       R_j(1,1) = -1;
       R_j(2,2) = 1;
       R_j(3,3) = 1;
    case 2
       R_j(1,1) = 1;
       R_j(2,2) = -1;
       R_j(3,3) = 1;
    case 3
       R_j(1,1) = 1;
       R_j(2,2) = 1;
       R_j(3,3) = -1;
end
