 % Preform a newton method based on a parameter and base point
% 
% Notes: 
% 
%% Clean MatLab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearAll
%% Map and Trajectory Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard Map Parameter alpha
alpha = pi/2;
%%%% Period 1
% initialP = [pi+.5+.5,0];
%%%% Period 6
initialP = [1.85, 0.565]; % large
% initialP = [1.3, -.75]; % small
%%%% Period 24
% initialP = [2.25, 1];
%%%% Period 3 wrapped around the edge (doesn't work yet)
% initialP = [pi, 3]; 
% alpha = pi/4;
% initialP = [pi,1];
%% Set up parameter variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Banach Space Parameter
nu = 1.1;
% Sobolev space H^m
sobolevMax = 10;
% Initial number of modes, and number of modes out of Newton
initialModes = 10;
currentModes = 50; 
modeStep = 20;
maxModes = 200;
% Number of points for plotting and computing initial parameterization
numPoints = 1000;
% Do we need to replace rho with 1-rho
rhoFlip = 1;
% Error limit on rho
rhoLimit = 13;
% Error limit in Newton
errorLimit = 6;
errorConjLimit = 6;
% Largest value allowed to change rho 
maxRhoStep = 10^-3;
% Do we need to subtract or add to rho 
delRhoFlip = 0; % Sign flip on delRho
% How small should be go on the steps 
maxPower = 20;
% Error limit in truncation
tailLimit = 16;
% % Attempt limit and related variables
maxAttempts = 50;
attemptNumber = 1;
successNumber = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% For general use, nothing below here needs to be changed.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up loop variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set total time to zero
totalTime = 0;
%% Ask for a file name for the diary output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = input('Name folder for output? (can be root/folder, etc.) ','s');
mkdir(fname);
diary(strcat('./',fname,'/output.txt'))
successFile = fopen(strcat('./', fname, '/successes.txt'), 'w');
%% Header information for the file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parameters: \n')
fprintf(['P =', repmat(' ',1,27),'[ %f, %f ]\n'], initialP(1),initialP(2))
fprintf(['alpha =', repmat(' ',1,23), '%.15g', repmat(' ',1,10),...
    'probably acos(something)\n'], alpha)
fprintf('nu =                          %f\n', nu)
fprintf('initial number of modes       %d\n', initialModes)
fprintf('number of modes after newton  %d\n', currentModes)
fprintf('number of points              %d\n', numPoints)
fprintf('maximum rho step         %.15g\n', maxRhoStep)
fprintf('minimum rho step         %.15g\n', maxRhoStep*10^-maxPower)
fprintf('max tail size            %.15g\n', 10^-tailLimit)
fprintf('max error                %.15g\n', 10^-errorLimit)
fprintf('rho difference limit     %.15g \n', 10^-rhoLimit)
fprintf('rho flip                 %d\n', rhoFlip)
fprintf('direction of rho change  %d\n', delRhoFlip)
fprintf(successFile, '\n'); 
fprintf(successFile, 'Parameters: \n');
fprintf(successFile, ['P =', repmat(' ',1,22),'[ %f, %f ]\n'], ...
    initialP(1),initialP(2));
fprintf(successFile, ['alpha =', repmat(' ',1,18), '%.15g', ...
    repmat(' ',1,20), 'probably acos(something)\n'], alpha);
fprintf(successFile, 'nu =                     %.15g\n', nu);
fprintf(successFile, 'sobolev space max m =    %d\n', sobolevMax);
fprintf(successFile, 'initial number of modes  %d\n', initialModes);
fprintf(successFile, 'Mode step used in newton %d\n', modeStep);
fprintf(successFile, 'number of points         %d\n', numPoints);
fprintf(successFile, 'maximum rho step         %.15g\n', maxRhoStep);
fprintf(successFile, 'minimum rho step         %.15g\n', ...
    maxRhoStep*10^-maxPower);
fprintf(successFile, 'max tail size            %.15g\n', 10^-tailLimit);
fprintf(successFile, 'max error                %.15g\n', 10^-errorLimit);
fprintf(successFile, 'rho difference limit     %.15g \n', 10^-rhoLimit);
fprintf(successFile, 'rho flip                 %d\n', rhoFlip);
fprintf(successFile, 'direction of rho change  %d\n', delRhoFlip);
fprintf(successFile, '\n');
%% Compute and plot trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Computing the phasespace and trajectory... \n    ')
tic
trajectory = pointTrajectory(initialP, alpha, numPoints);
toc
totalTime = totalTime + toc;
fprintf('\n')
%% Plot trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figTrajectory = figure;
hold on
set(figTrajectory, 'Units', 'Normalized', 'OuterPosition', ...
    [.05 .5 .25 .4]);
set(figTrajectory, 'Name', 'Trajectory and Initial Parameterization');
phasespacePlot(alpha);
plot(trajectory(1,:),trajectory(2,:),'.r');
axis([0, 2*pi, -4, 4]);
%% Ask the number of tori and separate the orbits %%%%%%%%%%%%%%%%%%%%%%%%%
K = input('How many periodic tori? ');
fprintf('\nSeparate the trajecotory... \n    ')
tic
tori = trajectorySeparator(trajectory, K);
toc
totalTime = totalTime + toc;
fprintf('\n')
%% Find periodic points of Tori %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Finding periodic points ... \n    ')
tic
perPts = perPointFinder(tori, K, alpha);
toc
totalTime = totalTime + toc;
fprintf('\n')
%% Plot periodic points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(perPts(1,:),perPts(2,:),'.k')
%% Find rotation number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Compute the rotation number... \n')
tic
rhoMultiplier = 1;
diffRho = 1;

numPointsRho = 0;
while diffRho > 10^-rhoLimit && numPointsRho < 10^8
    % Start with 1000 points per tori and add 1000 each iteration
    numPointsRho = K*1000*rhoMultiplier; 
    fprintf('    Number of points per tori         %d\n    ', 1000*rhoMultiplier)
    tic
    trajectory2 = pointTrajectory(initialP, alpha, numPointsRho);
    tori2 = trajectorySeparator(trajectory2, K);
    toc
    rho(rhoMultiplier) = weightedBirkoffRotationNumber(tori2, K, perPts);
    fprintf('    Rho approximates as %.15g. \n    ',rho(rhoMultiplier))
    toc
    if rhoMultiplier > 1
        diffRho = abs(rho(rhoMultiplier) - rho(rhoMultiplier-1));
        fprintf('    Difference between last rho and this rho: %.15d\n'...
            , abs(rho(rhoMultiplier)-rho(rhoMultiplier-1)))
        fprintf('    Machine epsilon is %.15d. \n\n', eps(rho(rhoMultiplier)))
    else % i <= 1
        fprintf('\n')
    end % end if
    rhoMultiplier = rhoMultiplier + 1;
end % end while loop
if numPointsRho >= 10^8
    fprintf('Rho accuracy questionable.\n')
end % end if
fprintf(successFile, 'The rotation number is approximately %.15g.\n    '...
    , rho(end));
fprintf('The rotation number is approximately %.15g.\n    ', rho(end))
rho = rho(end);
toc
totalTime = totalTime + toc;
fprintf('\n')
%% Initialize parameterization and pad with zeros %%%%%%%%%%%%%%%%%%%%%%%%%
param = cell(2,K);
fprintf('Compute the initial Fourier modes with %d initial modes... \n    ', initialModes)
tic
for k = 1:K
   [param{1,k}, param{2,k}] = fourierParam(tori{k}, initialModes, K, rho, perPts(:,k)); 
   param{1,k} = truncate(param{1,k}, currentModes);
   param{2,k} = truncate(param{2,k}, currentModes);
end % end for loop
toc
totalTime = totalTime + toc;
fprintf('\n')
%% Plot initial parameterization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figTrajectory);
fourierCellPlot(param, 'm');
figParameterization = figure;
fourierCellPlot(param, 'm');
%% Flip rho if needed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rhoFlip == 1 % rhoFlip(1) == '1'
    rho = 1 - rho;
    fprintf('    Flipped rho is %.15g.\n\n', rho)
else % rhoFlip == 0
    fprintf('\n')
end % end if
%% Measure the sequence space error of the initial parameterization %%%%%%%
fprintf('Compute initial parameterization and approx sine and cosine error with nu = %1.1f ... \n    ', nu)
tic
phase = zeros(2,K);
beta = 0;
scalars = zeros(2,K);
% Set up parameters and phase conditions
for k = 1:K
    param{3,k} = fourierSineApprox(param{1,k},50);
    param{4,k} = fourierCosineApprox(param{1,k},50);
    phase(1,k) = evaluate(param{1,k}, 0);
    phase(2,k) = evaluate(param{2,k}, 0);
end % end for loop
% Compute initial defect
defectError = normPhi(beta, scalars, param, alpha, rho, phase, nu);
fprintf('Error = %.15g \n     ', defectError)
toc
totalTime = totalTime + toc;
fprintf('\n')
%% Newton-like operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Carry out newton-like method...\n')
tic
newtonIter = 0;
while newtonIter < 10 || defectError >= 10^-errorLimit
   [betaNew, scalarsNew, paramNew] = newtonStep...
       (beta, scalars, param, alpha, rho, phase);
    newtonIter = newtonIter + 1;
    defectErrorNew = normPhi(betaNew, scalarsNew, paramNew, alpha, rho, phase, nu);
    if (defectErrorNew < defectError || newtonIter <= 2) && defectErrorNew < 50
        defectError = defectErrorNew;
        conjError = conjugacyError(rho, param, alpha, 1000);
        beta = betaNew;
        scalars = scalarsNew;
        param = paramNew;
        fprintf('    Iteration %d, error = %.15g \n', newtonIter, defectError)
        fprintf('       conjugacy error = %.15g \n     ', conjError)
        toc
    else % Reject new value
        fprintf('    Iteration %d, error exceeded limit \n     ',newtonIter)
        toc
        break
    end
end % while loop
fprintf('\n')
%% Compute sobolev norm and tail value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sobolevGrid = zeros(sobolevMax,1);
for m = 1:sobolevMax
    maxSobolevNorm = max(sobolevNorm(param,m));
    fprintf(successFile, 'max Sobolev m = %d Norm is %.15g \n', m, maxSobolevNorm);
    fprintf('max Sobolev m = %d Norm is %.15g \n', m, maxSobolevNorm)
    sobolevGrid(m, 1) = maxSobolevNorm;
end % for loop
endValue = maxEndVal(param);
fprintf(successFile, 'last coefficient absolute value is %.15g \n', endValue);
fprintf('last coefficient absolute value is %.15g \n', endValue)
totalTime = totalTime + toc;
fprintf('Total time elapsed %f seconds. \n', totalTime)
%% Plot newton output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figTrajectory);
fourierCellPlot(param, 'k')
%% Plot newton output on new figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figParameterization)
fourierCellPlot(param, 'k')
set(figParameterization, 'Units', 'Normalized', 'OuterPosition', [.25 .4 .25 .4]);
set(figParameterization, 'Name', 'Parameterizations')
axis([0,2*pi,-pi,pi])
%% Plot the log of the coefficents %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figCoefficients = figure;
set(figCoefficients, 'Units', 'Normalized', 'OuterPosition', [.47 .5 .25 .4]);
set(figCoefficients, 'Name', 'Log Plot of Coefficients');
% hold on
logCoeffPlot(param)
% logCoeffPlot(param{2},'k--')
% legend('x-param', 'y-param')
%% Plot the log of sobolev norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figSobolev = figure;
set(figSobolev, 'Units', 'Normalized', 'OuterPosition', [.67 .4 .25 .4]);
set(figSobolev, 'Name', 'Log Plot of Sobolev Norms');
semilogy(sobolevGrid,'b');
xlabel('Sobolev m value');
ylabel('Max norm')
%% Analytic continuation time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_original = rho;
paramOld = param;
%% Confirem contuation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input('\nContinue? (Ctrl-C to stop): ', 's');
currentPower = 0;
reRun = 0;
previousSuccess = param;
%% Main loop, shrinking rho size on unsuccessful attempts %%%%%%%%%%%%%%%%%
while currentPower <= maxPower && attemptNumber <= maxAttempts
    del_rho = ((-1)^delRhoFlip)*maxRhoStep*2^-currentPower;
    fprintf('\nUsing rho step of %.15g \n', del_rho)
    loopError = 0;
    % Iteration loop, attempting constant rho step size %%%%%%%%%%%%%%%%%%%
    while loopError < 10^-errorLimit && attemptNumber <= maxAttempts
        % Make rho smaller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rho = rho - del_rho;
        % Set up newton-like method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if reRun == 0
            fprintf('New rho =                               %.15f\n\n', rho)
        else % reRun == 1
            fprintf('Rho =                                   %.15f\n\n', rho)
        end % end if
        % Set up newton-like method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if reRun == 0
            paramOld = param;
        end
        fprintf('Start the new newton-like method using %d modes and nu = %1.1f ... \n', currentModes, nu)
        tic
        beta = 10^-K;
        for j = 1:K
            phase(1,j) = evaluate(param{1,j}, 0);
            phase(2,j) = evaluate(param{2,j}, 0);
        end
        loopError = normPhi(beta, scalars, param, alpha, rho, phase, nu);
        fprintf('    Iteration 0, error = %.15g \n     ', loopError)
        toc
        totalTime = totalTime + toc;

        % Carry out newton-like operator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        newtonIter = 0;
        while newtonIter < 10 || loopError >= 10^-errorLimit
           [betaNew, scalarsNew, paramNew] = newtonStep...
               (beta, scalars, param, alpha, rho, phase);
            newtonIter = newtonIter + 1;
           defectErrorNew = normPhi(betaNew, scalarsNew, paramNew, alpha, rho, phase, nu);
            if (defectErrorNew < loopError || newtonIter <= 2) && defectErrorNew < 50
                loopError = defectErrorNew;
                conjError = conjugacyError(rho, param, alpha, 1000);
                beta = betaNew;
                scalars = scalarsNew;
                param = paramNew;
                fprintf('    Iteration %d, error = %.15g \n', newtonIter, loopError)
                fprintf('       conjugacy error = %.15g \n     ', conjError)
                toc
            else % Reject new value
                fprintf('    Iteration %d, error exceeded limit \n     ',newtonIter)
                toc
                break
            end
        end % while loop
        fprintf('\n')        
        % Compute output metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        endValue = maxEndVal(param);
        fprintf('The end coefficients absolute value is: %.15g \n', endValue)
        totalTime = totalTime + toc;
        fprintf('    Total time elapsed %f seconds. \n', totalTime)
        % Plot and save initial coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(figCoefficients)
        logCoeffPlot(param)
        % Decide next move: continue, rerun with more modes, %%%%%%%%%%%%%%
        % rerun with smaller rhostep, stop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Accept new parameterization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if endValue < 10^-tailLimit && loopError < 10^-errorLimit && conjError < 10^-errorConjLimit
            figure(figParameterization)
            fourierCellPlot(param)
            previousSuccess = param;
            % Record success            
            successNumber = successNumber + 1;
            fprintf(successFile, 'Success number %d \n', successNumber);
            fprintf(successFile, 'Attempt number %d, successful!\n', attemptNumber);
            fprintf(successFile, '    number of modes    %d\n', currentModes);
            fprintf(successFile, '    error              %.15g\n', loopError);
            fprintf(successFile, '    tail value         %.15g\n', endValue);
            fprintf(successFile, '    rho                %f\n', rho);
            % Compute Sobolev norms
            tempSobVect = zeros(sobolevMax, 1);
            for m = 1:sobolevMax
                maxSobolevNorm = max(sobolevNorm(param,m));
                fprintf('The max Sobolev m = %d Norm is:                %.15g \n', m, maxSobolevNorm)
                fprintf(successFile, '    sobolev norm m = %d   %.15g\n', m, maxSobolevNorm);
                tempSobVect(m, 1) = maxSobolevNorm;
            end
            sobolevGrid = [sobolevGrid, tempSobVect];
            fprintf(successFile, '\n');
            fprintf('Attempt number %d, successful!\n', attemptNumber)
            fprintf('Success number %d!\n', successNumber)
            attemptNumber = attemptNumber + 1;
            % Truncate back to original number of modes and exit loop
            currentModes = initialModes;
            for k = 1:K
                param{1,k} = truncate(param{1,k},currentModes);
                param{2,k} = truncate(param{2,k},currentModes);
                param{3,k} = truncate(param{3,k},currentModes);
                param{4,k} = truncate(param{4,k},currentModes);
            end
            % Reset stepsize of del_rho on a success
            currentPower = 0;
            loopError = 1;
            reRun = 0;
        % Rerun with more modes or move on if too many modes %%%%%%%%%%%%%%
        elseif endValue > 10^-tailLimit
            rho = rho + del_rho;
            currentModes = currentModes + modeStep;
            for k = 1:K
                param{1,k} = truncate(param{1,k},currentModes);
                param{2,k} = truncate(param{2,k},currentModes);
                param{3,k} = truncate(param{3,k},currentModes);
                param{4,k} = truncate(param{4,k},currentModes);
            end
            if currentModes <= maxModes
                loopError = 0;
                fprintf('\nMore modes required, rerunning Newton with %d modes.\n', currentModes)
                reRun = 1;
            else
                if currentPower < maxPower 
                    fprintf('\nMax modes exceeded, trying again with smaller rho step.\n')
                else % i >= maxPower
                    fprintf('\nMax modes exceeded, and already at smallest rho step.\n')
                end % if
                currentModes = initialModes;
                param = paramOld;
                for k = 1:K
                    param{1,k} = truncate(param{1,k},currentModes);
                    param{2,k} = truncate(param{2,k},currentModes);
                    param{3,k} = truncate(param{3,k},currentModes);
                    param{4,k} = truncate(param{4,k},currentModes);
                end
                reRun = 0;
                currentPower = currentPower + 1;
                fprintf('\nAttempt number %d, failed!\n', attemptNumber)
                attemptNumber = attemptNumber + 1;
                loopError = errorLimit;
            end
        % Try smaller stepsize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            currentModes = initialModes;
            param = paramOld;
            for k = 1:K
                param{1,k} = truncate(param{1,k},currentModes);
                param{2,k} = truncate(param{2,k},currentModes);
                param{3,k} = truncate(param{3,k},currentModes);
                param{4,k} = truncate(param{4,k},currentModes);
            end
            reRun = 0;
            rho = rho + del_rho;
            currentPower = currentPower + 1;
            fprintf('\nAttempt number %d, failed!\n', attemptNumber)
            attemptNumber = attemptNumber + 1;
        end
    end
end
% Plot and save initial coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figCoefficients)
logCoeffPlot(previousSuccess)
%% Plot the sobolev norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figSobolev)
if size(sobolevGrid,2) == 1
    sobolevGrid = [sobolevGrid, zeros(sobolevMax,1)];
end
mesh(log(sobolevGrid)','FaceColor','flat', 'FaceAlpha', '0.75')
axis tight
xlabel('Sobolev m value')
ylabel('parameterization number')
zlabel('norm')
%% Save figures and close stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(figTrajectory, strcat('./',fname,'/phaseSpace.png'));
saveas(figParameterization, strcat('./',fname,'/ac_final_parameterization.png'));
saveas(figCoefficients,strcat('./',fname,'/ac_final_logCoeffs.png'));
saveas(figSobolev,strcat('./',fname,'/sobolevNorms.png'));
fclose(successFile);
diary off
fclose all;