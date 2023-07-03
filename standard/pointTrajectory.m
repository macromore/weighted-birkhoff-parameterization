function orbit = pointTrajectory(P, alpha, numPoints)
%Given an initial point P and a parameter create and plot a trajectory of
%the point. 
% 
% Define variables:
% P             - Base point (input)
% alpha         - standard map parameter (input)
% numPoints     - length of trajectory (input)
% s             - dummy variable to set up P 
% thisPoint     - temp point variable
% orbit         - trajectory (output)
% 
% Dependencies: 
% standardMap.m
    if nargin == 1
        alpha = pi/2;
        numPoints = 250;
    elseif nargin == 2 
        numPoints = 250;
    end % end if
    s = size(P,1);
    if s(1) == 1
        P = P';
    end % end if
    thisPoint = [P; 0];
    orbit = zeros(2, numPoints);
    orbit(1,1) = thisPoint(1);
    orbit(2,1) = thisPoint(2);
    for k = 2:numPoints
       thisPoint = standardMap(thisPoint, alpha);
       orbit(:,k) = thisPoint;
    end % end for loop
end % end pointTrajectory