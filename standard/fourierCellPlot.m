function fourierCellPlot(param, color)
% Plot multiple pairs of Fourier parameterizations which are stored in a
% cell array
%
% Define variables:
% param         - cell array containing pairs of Fourier series
% holdState     - stores the value of plot hold
% 
% Dependencies: 
% pairPlot.m
% Fourier.m
    % Check current hold status and turn it on
    holdState = ishold;
    if holdState == 0
        hold on
    end % end if
    % Plot the pairs of Fourier series
    if nargin == 2
        for i = 1:size(param,2)
            pairPlot(param{1,i},param{2,i}, color);            
        end % end for loop
    else % nargin != 2
        for i = 1:size(param,2)
            pairPlot(param{1,i},param{2,i});            
        end % end for loop
    end % end if
    %Turn hold off if it was off to begin with
    if holdState == 0
        hold off
    end % end if
end % end