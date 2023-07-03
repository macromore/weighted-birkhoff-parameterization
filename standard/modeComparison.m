%% Clean MatLab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearAll
%% Map and Trajectory Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard Map Parameter alpha
%%%% Period 1
% initialP = [pi+.5,0];
%%%% Period 6
% initialP = [1.85, 0.565]; % large
% initialP = [1.3, -.75]; % small
%%%% Period 24
% initialP = [2.25, 1];
% initialP = [5, .5];
%%%% Period 3 wrapped around the edge (doesn't work yet)
% initialP = [pi, 3]; 
% alpha = pi/4;
% initialP = [pi,1];
alpha = pi/2;
% initialP = [4.85, 0.5];
% initialP = [4.8, 0.5];
initialP = [4.8155, 0.5];
%% Set up parameter variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Banach Space Parameter
nu = 1.1;
% Sobolev space H^m
sobolevMax = 10;
% Conjugacy iteration max
conjMax = 100;
% Initial number of modes, and number of modes out of Newton
initialModes = 25;
modeStep = 25;
maxModes = 5;
% Number of points for plotting and computing initial parameterization
numPoints = 24000;
% Do we need to replace rho with 1-rho
rhoFlip = 1;
% Error limit on rho
rhoLimit = 15;
% Error limit in Newtons
errorLimit = 10;
%% Plot color %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color = 'm';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% For general use, nothing below here needs to be changed.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ask for a file name for the diary output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = input('Name folder for output? (can be root/folder, etc.) ','s');
mkdir(fname);
diary(strcat('./',fname,'/output.txt'))
tableFile = fopen(strcat('./', fname, '/table.txt'), 'w');
%% Header information for the file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parameters: \n')
fprintf(['P =', repmat(' ',1,27),'[ %f, %f ]\n'], initialP(1),initialP(2))
fprintf(['alpha =', repmat(' ',1,23), '%.15g', repmat(' ',1,10),...
    'probably acos(something)\n'], alpha)
fprintf('nu =                      %f\n', nu)
fprintf('sobolev max space m =    %d\n', sobolevMax)
fprintf('initial number of modes  %d\n', initialModes)
fprintf('Mode step                %d\n', modeStep)
fprintf('points per circle        %d\n', numPoints)
fprintf('max error                %.15g\n', 10^-errorLimit)
fprintf('rho difference limit     %.15g \n', 10^-rhoLimit)
fprintf('rho flip                 %d\n', rhoFlip)
fprintf(tableFile, '\\begin{tabular}{c||c|c|c|c}\n Numer of Modes & Newton Iterations & Tail Size & $\\|\\Psi\\|$ Error & Conjugacy Error\\\\ \n \\hline \n');
%% Compute and plot trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Computing the phasespace and trajectory... \n  ')
tic
trajectory = pointTrajectory(initialP, alpha, numPoints);
toc
fprintf('\n')
%% Plot trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figTrajectory = figure(1);
set(figTrajectory, 'Units', 'Normalized', 'OuterPosition', ...
    [.05 .5 .25 .4]);
set(figTrajectory, 'Name', 'Trajectory and Initial Parameterization');
hold on
phasespacePlot(alpha);
plot(trajectory(1,:),trajectory(2,:),'.r');
axis([0, 2*pi, -4, 4]);
%% Ask the number of tori and separate the orbits %%%%%%%%%%%%%%%%%%%%%%%%%
K = input('How many periodic tori? ');
fprintf('\nSeparate the trajecotory... \n  ')
tori = trajectorySeparator(trajectory, K);
toc
fprintf('\n')
%% Find periodic points of Tori %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Finding periodic points ... \n  ')
perPts = perPointFinder(tori, K, alpha);
toc
fprintf('\n')
figure(1)
plot(perPts(1,:),perPts(2,:),'.k')
%% Find rotation number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Compute the rotation number... \n')
rhoMultiplier = 1;
diffRho = 1;
numPointsRho = 0;
while diffRho > 10^-rhoLimit && numPointsRho < 10^8
    % Start with 1000 points per tori and add 1000 each iteration
    numPointsRho = K*1000*rhoMultiplier; 
    fprintf('Number of points per tori         %d\n', 1000*rhoMultiplier)
    trajectory2 = pointTrajectory(initialP, alpha, numPointsRho);
    tori2 = trajectorySeparator(trajectory2, K);
    rho(rhoMultiplier) = weightedBirkoffRotationNumber(tori2, K, perPts);
    fprintf('    Rho approximates as %.15g. \n',rho(rhoMultiplier))
    if rhoMultiplier > 1
        diffRho = abs(rho(rhoMultiplier) - rho(rhoMultiplier-1));
        fprintf('  Delta rho:        %e \n'...
            , abs(rho(rhoMultiplier)-rho(rhoMultiplier-1)))
    end % end if
    fprintf('  Machine epsilon:  %e \n    ', eps(rho(rhoMultiplier)))
    toc
    rhoMultiplier = rhoMultiplier + 1;
end % end while loop
if numPointsRho >= 10^8
    fprintf('Rho accuracy questionable.\n')
end % end if
fprintf('The rotation number is approximately %.15g.\n  ', rho(end))
rho = rho(end);
toc
fprintf('\n')
%% Initialize parameterization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param = cell(2,K);
fprintf('Compute the initial Fourier modes with %d initial modes... \n  ', initialModes)
for k = 1:K
   [param{1,k}, param{2,k}] = fourierParam(tori{k}, initialModes, K, rho, perPts(:,k)); 
end % end for loop
toc
fprintf('\n')
%% Plot initial parameterization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figTrajectory);
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
defectError = normPhi(beta, scalars, param, alpha, rho, phase, nu);
conjError = conjugacyError(rho, param, alpha, conjMax);
fprintf('Error           = %e \n     ', defectError)
fprintf('Conjugacy Error = %e \n', conjError)
toc
fprintf('\n')
paramOrig = param;
sobolevGrid = zeros(sobolevMax,1);
figCoefficients = figure(2);
logCoeffPlot(param)
set(figCoefficients, 'Units', 'Normalized', 'OuterPosition', [.47 .5 .25 .4]);
set(figCoefficients, 'Name', 'Log Plot of Coefficients');
%% Set Up plot for sobolev norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figSobolev = figure(3);
set(figSobolev, 'Units', 'Normalized', 'OuterPosition', [.67 .4 .25 .4]);
set(figSobolev, 'Name', 'Log Plot of Sobolev Norms');
%% Loop to add modes
for currentModes = initialModes:modeStep:maxModes
    fprintf('\nNumer of modes : %d\n', currentModes)
    % Truncate to current modes
    for k = 1:K
        for j = 1:4
            param{j,k} = truncate(paramOrig{j,k}, currentModes);
        end
    end
    % Newton-like operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('  Carry out newton-like method...\n')
    newtonIter = 0;
    while newtonIter < 10 || defectError >= 10^-errorLimit
        [betaNew, scalarsNew, paramNew] = newtonStep...
            (beta, scalars, param, alpha, rho, phase);
        defectErrorNew = normPhi(betaNew, scalarsNew, paramNew, alpha, rho, phase, nu);
        if defectErrorNew < defectError || newtonIter < 1
            newtonIter = newtonIter + 1;
            defectError = defectErrorNew;
            beta = betaNew;
            scalars = scalarsNew;
            param = paramNew;
            conjError = conjugacyError(rho, param, alpha, conjMax);
                fprintf('    Iteration %d, error = %e \n', newtonIter, defectError)
                fprintf('       conjugacy error = %e \n', conjError)
        else % Reject new value
            fprintf('    Iteration %d, error limit exceeded \n',newtonIter + 1)
            if newtonIter == 0
                fprint('  Newton Method Failed.\n\n')
            end
            break
        end
    end % while loop
    % Compute sobolev norm and tail value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sobolevTemp = zeros(sobolevMax,1);
    fprintf('  Sobolev Norms:\n');
    for m = 1:sobolevMax
        maxSobolevNorm = max(sobolevNorm(param,m));
        fprintf('    k = %d norm is %e \n', m, maxSobolevNorm)
        sobolevTemp(m) = maxSobolevNorm;
    end % for loop
    sobolevGrid = [sobolevGrid sobolevTemp];
    endValue = maxEndVal(param);
    fprintf('  Tail size: %e \n', endValue)
    fprintf(tableFile, ' %d & %d & %e & %e & %e \\\\ \n', currentModes, newtonIter, endValue, defectError, conjError);
    % Plot newton output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(figTrajectory);
    fourierCellPlot(param, 'k')
    % Plot the log of the coefficents %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(figCoefficients);
    logCoeffPlot(param)
    axis tight
    % Plot sobolev norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if currentModes > initialModes
        figure(figSobolev)
        mesh(1:sobolevMax, initialModes:modeStep:currentModes, log(sobolevGrid(:,2:end))','FaceColor','flat', 'FaceAlpha', '0.75')
        axis tight        
        xlabel('Sobolev k value')
        ylabel('Numer of Modes')
        zlabel('Log of norm')
    end
    toc
end
%% Save figures and close stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(figTrajectory, strcat('./',fname,'/phaseSpace.png'));
saveas(figCoefficients,strcat('./',fname,'/logCoeffs.png'));
saveas(figSobolev,strcat('./',fname,'/sobolevNorms.png'));
fprintf(tableFile, '\\end{tabular} %% Matlab generated table');
fclose(tableFile);
diary off;
fclose all;