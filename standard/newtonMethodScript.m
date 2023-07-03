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
% initialP = [pi+.5,0];
%%%% Period 6
initialP = [1.85, 0.565]; % large
% initialP = [1.3, -.75]; % small
%%%% Period 24
% initialP = [2.25, 1];
% initialP = [4.8, .5];
% initialP = [4.8155, 0.5];
%%%% Period 3 wrapped around the edge (doesn't work yet)
% initialP = [pi, 3]; 
% alpha = pi/4;
% initialP = [pi,1];
%% Set up parameter variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Banach Space Parameter
nu = 1.1;
% Sobolev space H^m
sobolevMax = 10;
% Conjugacy check max
conjMax = 100;
% Initial number of modes, and number of modes out of Newton
initialModes = 20;
newtonModes = 50; 
% Number of points for plotting and computing initial parameterization
numPoints = 24000;
% Do we need to replace rho with 1-rho
rhoFlip = 1;
% Error limit on rho
rhoLimit = 15;
% Error limit in Newton
errorLimit = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% For general use, nothing below here needs to be changed.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ask for a file name for the diary output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = input('Name folder for output? (can be root/folder, etc.) ','s');
mkdir(fname);
diary(strcat('./',fname,'/output.txt'))
%% Header information for the file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parameters: \n')
fprintf(['P =', repmat(' ',1,27),'[ %f, %f ]\n'], initialP(1),initialP(2))
fprintf(['alpha =', repmat(' ',1,23), '%.15g', repmat(' ',1,10),...
    'probably acos(something)\n'], alpha)
fprintf('nu =                          %f\n', nu)
fprintf('initial number of modes       %d\n', initialModes)
fprintf('number of modes after newton  %d\n', newtonModes)
fprintf('number of points              %d\n\n', numPoints)
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
fprintf('\nSeparate the trajecotory... \n    ')
tori = trajectorySeparator(trajectory, K);
toc
fprintf('\n')
%% Find periodic points of Tori %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Finding periodic points ... \n  ')
perPts = perPointFinder(tori, K, alpha);
toc
fprintf('\n')
plot(perPts(1,:),perPts(2,:),'.k')
%% Find rotation number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Compute the rotation number... \n')
rhoMultiplier = 10;
diffRho = 1;
numPointsRho = 0;
while diffRho > 10^-rhoLimit && numPointsRho < 10^8
    % Start with 1000 points per tori and add 1000 each iteration
    numPointsRho = K*1000*rhoMultiplier; 
    fprintf('Number of points per tori         %d\n    ', 1000*rhoMultiplier)
    trajectory2 = pointTrajectory(initialP, alpha, numPointsRho);
    tori2 = trajectorySeparator(trajectory2, K);
    rho(rhoMultiplier) = weightedBirkoffRotationNumber(tori2, K, perPts);
    fprintf('  Rho approximates as %.15g. \n',rho(rhoMultiplier))
    if rhoMultiplier > 10
        diffRho = abs(rho(rhoMultiplier) - rho(rhoMultiplier-1));
        fprintf('  Delta rho:          %e \n'...
            , abs(rho(rhoMultiplier)-rho(rhoMultiplier-1)))
    end % end if
    fprintf('  Machine epsilon is  %e \n  ', eps(rho(rhoMultiplier)))
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
%% Initialize parameterization and pad with zeros %%%%%%%%%%%%%%%%%%%%%%%%%
param = cell(2,K);
fprintf('Compute the initial Fourier modes with %d initial modes... \n  ', initialModes)
for k = 1:K
   [param{1,k}, param{2,k}] = fourierParam(tori{k}, initialModes, K, rho, perPts(:,k)); 
   param{1,k} = truncate(param{1,k}, newtonModes);
   param{2,k} = truncate(param{2,k}, newtonModes);
end % end for loop
toc
fprintf('\n')
%% Plot initial parameterization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figTrajectory);
fourierCellPlot(param, 'm');
figParameterization = figure(2);
fourierCellPlot(param, 'm');
paramOrig = param;
%% Flip rho if needed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rhoFlip == 1 % rhoFlip(1) == '1'
    rho = 1 - rho;
    fprintf('Flipped rho is %.15g.\n\n', rho)
else % rhoFlip == 0
    fprintf('\n')
end % end if
%% Measure the sequence space error of the initial parameterization %%%%%%%
fprintf('Compute initial parameterization and approx sine and cosine error with nu = %1.1f ... \n', nu)
phase = zeros(2,K);
beta = 0;
scalars = zeros(2,K);
% Set up parameters and phase conditions
for k = 1:K
    param{3,k} = fourierSineApprox(param{1,k},20);
    param{4,k} = fourierCosineApprox(param{1,k},20);
    phase(1,k) = evaluate(param{1,k}, 0);
    phase(2,k) = evaluate(param{2,k}, 0);
end % end for loop
% Compute initial defect
defectError = normPhi(beta, scalars, param, alpha, rho, phase, nu);
conjError = conjugacyError(rho, param, alpha, conjMax);
fprintf('Error           = %e \n  ', defectError)
fprintf('Conjugacy Error = %e \n  ', conjError)
toc
fprintf('\n')
%% Newton-like operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Carry out newton-like method...\n')
newtonIter = 0;
while newtonIter < 10 || defectError >= 10^-errorLimit
   [betaNew, scalarsNew, paramNew] = newtonStep...
       (beta, scalars, param, alpha, rho, phase);
   defectErrorNew = normPhi(betaNew, scalarsNew, paramNew, alpha, rho, phase, nu);
   if defectErrorNew < defectError % || newtonIter <= 2) && defectErrorNew < 50
        newtonIter = newtonIter + 1;
        defectError = defectErrorNew;
        beta = betaNew;
        scalars = scalarsNew;
        param = paramNew;
        conjError = conjugacyError(rho, param, alpha, 1000);
        fprintf('    Iteration %d, error = %e \n', newtonIter, defectError)
%         fprintf('    Machine Epsilon:  %e \n', eps(defectError))
        fprintf('       conjugacy error = %e \n', conjError)
    else % Reject new value
            fprintf('    Iteration %d, error limit exceeded \n',newtonIter + 1)
        break
    end
end % while loop
fprintf('\n')
%% Compute sobolev norm and tail value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sobolevGrid = zeros(sobolevMax,1);
    fprintf('  Sobolev Norms:\n');
for m = 1:sobolevMax
    maxSobolevNorm = max(sobolevNorm(param,m));
    fprintf('    k = %d norm is %e \n', m, maxSobolevNorm)
    sobolevGrid(m, 1) = maxSobolevNorm;
end % for loop
endValue = maxEndVal(param);
fprintf('  Tail size: %e \n', endValue)
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
figCoefficients = figure(3);
set(figCoefficients, 'Units', 'Normalized', 'OuterPosition', [.47 .5 .25 .4]);
set(figCoefficients, 'Name', 'Log Plot of Coefficients');
hold on
logCoeffPlot(param)
axis tight
% logCoeffPlot(param{2},'k--')
% legend('x-param', 'y-param')
%% Plot the log of sobolev norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figSobolev = figure(4);
set(figSobolev, 'Units', 'Normalized', 'OuterPosition', [.67 .4 .25 .4]);
set(figSobolev, 'Name', 'Log Plot of Sobolev Norms');
semilogy(sobolevGrid,'b');
axis tight
xlabel('Sobolev m value');
ylabel('Max norm')
%% Save figures and close stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(figTrajectory, strcat('./',fname,'/phaseSpace.png'));
saveas(figParameterization, strcat('./',fname,'/parameterization.png'));
saveas(figCoefficients,strcat('./',fname,'/logCoeffs.png'));
saveas(figSobolev,strcat('./',fname,'/sobolevNorms.png'));
diary off;