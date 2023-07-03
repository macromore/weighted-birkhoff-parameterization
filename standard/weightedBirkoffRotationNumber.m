function rho = weightedBirkoffRotationNumber(orbit, K, perPts)
% Given a periodic trajectory, compute a rotation number using weighted
% averages
%
% Define variables:
% orbit         - initial data, the trajectory (input)
% K             - number of closed orbits (input)
% rho           - computed rotation number
% avex          - average of the x trajectory
% avey          - average of the y trajectory
% numPoints     - number of iterates to use
% args          - angle coordinates of points
% argDiff       - change in angle between two points
% totalWeight   - denominator for weight function
% weightValue   - numerator for weight function
% partialSum    - partial sum of weighted average
% 
% Dependencies: 
% N/A
    % set rho to zero as a trigger for recursion
    rho = 0;
    if K ~= 1
       rho = zeros(1,K);
       for i = 1:K
          rho(i) = weightedBirkoffRotationNumber(orbit{i}, 1, perPts(:,i));
       end % end for loop
       rho = median(rho)/K;
    end % end if
    % if we haven't computed a rotation number yet, compute rotation number
    if rho == 0
        if isa(orbit,'cell')
           orbit = orbit{1}; 
        end % end if
        % translate torus to origin
        orbit(1,:) = orbit(1,:) - perPts(1);
        orbit(2,:) = orbit(2,:) - perPts(2);
        args = atan2(orbit(2,:),orbit(1,:))/(2*pi);
        argDiff = mod(diff(args),1);
        % Compute the weighted sum
        weights = vectorWeight(length(orbit));
        rho = sum(weights(1,2:end) .* argDiff);
    end % end if
end % end weighted BirkhoffRotationNumber