function tplot = trajectoryPlot(P,alpha,numPoints,color)
%Given a point, a parameter, and a color, produce a plot of the trajectory.
% 
% Define variables:
% P             - base point (input)
% alpha         - standard map parameter (input)
% numPoints     - number of points to compute (input)
% color         - desired graph color (input)
% orbit         - orbit of the point 
% ty            - hold state 
% tplot         - plot of trajectory (output)
% 
% Dependencies: 
% pointTrajectory.m
    if nargin == 1
        alpha = pi/2;
        numPoints = 500;
        color = 'b';    
    elseif nargin == 2
        numPoints = 500;
        color = 'b';
    elseif nargin == 3
        color = 'b';
    end % end if
    orbit = pointTrajectory(P, alpha, numPoints);
    if size(findall(0,'type','figure')) ~= 0
        ty = ishold;
    else % only one figure
        ty = 0;
    end % end if
    if ty == 0
        figure
        hold on
    end % end if
    tplot = plot(orbit(1, :), orbit(2,:), strcat(color,'.'));
    axis([0, 2*pi, -pi, pi])
    if ty == 0
        hold off
    end % end if
end % end trajectoryPlot