function ppoints = perPointFinder(orbits, K, alpha)
% Find periodic points inside a trajectory
%
% Define variables:
% tori      - separated trajectory (input)
% K         - number of orbits (input)
% alpha     - standard map parameter (input)
% ppoints   - approximate periodic points (output)
% 
% Dependencies: 
% N/A
    ppoints = zeros(2,K);
    rho = 1;
    param = cell(4,K);
    beta = 0;
    derr = 1;
    err = 1;
    iterationCount = 1;
    phase = zeros(1,K);
    scalars = zeros(2,K);
    for k = 1:K
        % Use average of orbit as initial guess
        ppoints(:,k) = mean(orbits{k},2);
        % Convert to Fourier series
        param{1,k} = Fourier(ppoints(1,k));
        param{2,k} = Fourier(ppoints(2,k));
        param{3,k} = fourierSineApprox(param{1,k},50);
        param{4,k} = fourierCosineApprox(param{1,k},50);
%         phase(1,k) = evaluate(param{1,k}, 0);
%         phase(2,k) = evaluate(param{2,k}, 0);
    end % for loop
    plot(ppoints(1,:),ppoints(2,:),'.m')
    while  0 == 1 % derr > 10^-5 && iterationCount < 10
%         for k = 1:K
%             phase(1,k) = evaluate(param{1,k}, 0);
%         end % end for loop
        [beta, scalars, param] = ...
            newtonStep(beta, scalars, param, alpha, rho, phase);
        temp = [evaluate(param{1,1},0); evaluate(param{2,1},0)];
        for dummyVar = 1:K
            temp = standardMap(temp, alpha);
        end % end for loop
        err1 = norm([evaluate(param{1,1},0);evaluate(param{2,1},0)] - temp);
        derr = abs(err - err1);
        err = err1;
        iterationCount = iterationCount + 1;
    end % end while loop
    fprintf('Periodic point conjugacy error is %1.15f\n', err)
    for k = 1:K   
       ppoints(1,k) = evaluate(param{1,k},0);
       ppoints(2,k) = evaluate(param{2,k},0);   
    end
end % perPointFinder