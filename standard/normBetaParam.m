function val = normBetaParam(beta, param, nu)
% Compute the sequence space norm of the parameterization
%
% Define variables:
% beta      - beta from map (input)
% param     - parameterization (input)
% nu        - sequence space parameter (input)
% val       - norm value (output)
% 
% Dependencies: 
% Fourier.m
    val = abs(beta);
    for k = 1:size(param,2)
% Lazy coding should fix
       val = [val, norm(param{1,k}, nu), norm(param{2,k}, nu)]; 
    end % end for loop
    val = max(val);
end % end normBetaParam