function [betanext, scalarsnext, paramnext] = phi...
    (beta, scalars, param, alpha, rho, phase)
% the phi map
% 
% Define variables:
% beta          - beta from the map (input)
% scalars       - scalars for the map (intput)
% param         - cell array containing pairs of Fourier series (input)
% alpha         - standard map parameter (input)
% rho           - rotation number (input)
% phase         - vector of phase conditions (input)
% betanext      - beta (output)
% scalarsnext   - scalars (output)
% paramnext     - parameters (output)
% K             - number of modes
% 
% Dependencies: 
% Fourier.m
    K = size(param,2);
    betanext = evaluate(param{1,1}, 0) - phase(1,1);
    scalarsnext = zeros(2,K);
    paramnext = cell(4,K);
% Here is a possible issue we are locking the scalars down to the original phase
    for k = 1:K
        scalarsnext(1,k) = evaluate(param{3,k}, 0) - sin(phase(1,k));
        scalarsnext(2,k) = evaluate(param{4,k}, 0) - cos(phase(1,k));
    end % end for loop
    for k = 1:K-1
        paramnext{1,k} = param{1,k} + param{2,k} + alpha*param{3,k}...
            - param{1,k+1};
        paramnext{2,k} = param{2,k} + alpha*param{3,k} - param{2,k+1};
        paramnext{3,k} = diff(param{3,k}) - param{4,k}*diff(param{1,k})...
            + scalars(1,k)*param{3,k} - scalars(2,k)*param{4,k};
        paramnext{4,k} = diff(param{4,k}) + param{3,k}*diff(param{1,k})...
            + scalars(1,k)*param{4,k} + scalars(2,k)*param{3,k};
    end % end for loop
    paramnext{1,K} = param{1,K} + param{2,K} + alpha*param{3,K} - ...
        (1+beta)*rotation(param{1,1},rho*K);
    paramnext{2,K} = param{2,K} + alpha*param{3,K} - ...
        (1+beta)*rotation(param{2,1}, rho*K);
    paramnext{3,K} = diff(param{3,K}) - param{4,K}*diff(param{1,K})...
        + scalars(1,K)*param{3,K} - scalars(2,K)*param{4,K};
    paramnext{4,K} = diff(param{4,K}) + param{3,K}*diff(param{1,K})...
        + scalars(1,K)*param{4,K} + scalars(2,K)*param{3,K};
end % end phi