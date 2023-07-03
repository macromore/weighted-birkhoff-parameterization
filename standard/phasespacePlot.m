function phasespacePlot(alpha,numTrajectories, numPoints, color)
% plot a few trajectories
% 
% Define variables:
% P         - point to generate trajectory (input)
% alpha     - area preserving Henon map parameter (input)
% numPoints - length of trajectory to plot (input)
% color     - color of plot (input)
% 
% Dependencies: 
% trajectoryPlot.m
    %Generate and plot orbits
    if nargin == 0
        alpha = pi/2; 
        numTrajectories = 100;
        numPoints = 250;
        color = 'k';
    elseif nargin == 1
        numTrajectories = 100;
        numPoints = 250;
        color = 'k';
    elseif nargin == 2
        numPoints = 250;
        color = 'k';
    elseif nargin == 3
        color = 'k';
    end % end if
    P = linspace(-pi,pi, numTrajectories);
    if size(findall(0,'type','figure')) ~= 0
        ty = ishold;
    else % only one figure
        ty = 0;
    end % end if
    if ty == 0
        figure
        hold on
    end % end if
    for n = 1:numTrajectories
        trajectoryPlot([pi,P(n)], alpha, numPoints, color);
    end % end for loop
    axis([0, 2*pi, -pi, pi])
    if ty == 0
        hold off
    end % end if
end %end phasespacePlot