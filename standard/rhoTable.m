% Compute the same circle with successively more modes and report the
% results.
% 
% Notes: 
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
alpha = pi/4;
initialP = [pi,1];
% alpha = pi/2;
% initialP = [4.85, 0.5];
% initialP = [4.8, 0.5];
% initialP = [4.8155, 0.5];
%% Set up parameter variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Banach Space Parameter
nu = 1.1;
% Sobolev space H^m
sobolevMax = 10;
% Conjugacy iteration max
conjMax = 2;
% Initial number of modes, and number of modes out of Newton
initialModes = 25;
modeStep = 25;
maxModes = 350;
% Number of points for plotting and computing initial parameterization
numPoints = 24000;
% Do we need to replace rho with 1-rho
rhoFlip = 1;
% Error limit on rho
rhoLimit = 15;
% Error limit in Newtons
errorLimit = 10;

tailLimit = 15;
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
fprintf(['P =', repmat(' ',1,22),'[ %f, %f ]\n'], initialP(1),initialP(2))
fprintf(['alpha =', repmat(' ',1,18), '%.15g', repmat(' ',1,10),...
    'probably acos(something)\n'], alpha)
fprintf('nu =                     %.15g\n', nu)
fprintf('sobolev max space k      %d\n', sobolevMax)
fprintf('initial number of modes  %d\n', initialModes)
fprintf('Mode step                %d\n', modeStep)
fprintf('points per circle        %d\n', numPoints)
fprintf('max tail size            %.15g\n', 10^-tailLimit)
fprintf('max error                %.15g\n', 10^-errorLimit)
fprintf('Conjugacy iterations     %d\n', conjMax)
% fprintf('rho difference limit     %.15g \n', 10^-rhoLimit)
fprintf('rho flip                 %d\n\n', rhoFlip)
fprintf(tableFile, '\\begin{tabular}{c||l|c|c}\n Numer of Points &  $\\rho$ Approximation & Relative Error & Machine Epsilon \\\\ \n \\hline \n');
%% Compute trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Computing the trajectory... \n  ')
tic
trajectory = pointTrajectory(initialP, alpha, numPoints);
toc
fprintf('\n')
%% Plot trajectory and figure hold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figTrajectory = figure(1);
set(figTrajectory, 'Units', 'Normalized', 'OuterPosition', ...
    [.05 .5 .25 .4]);
set(figTrajectory, 'Name', 'Trajectory and Initial Parameterization');
hold on
phasespacePlot(alpha)
plot(trajectory(1,:),trajectory(2,:),'.g');
% axis([-2 2 -2 2]) % Do I need to do this?
%% Ask the number of tori and separate the orbits %%%%%%%%%%%%%%%%%%%%%%%%%
K = input('How many periodic tori? ');
fprintf('\nSeparate the trajecotory... \n    '    )
trajectory = pointTrajectory(initialP, alpha, numPoints*K);
tori = trajectorySeparator(trajectory, K);
toc
fprintf('\n')
%% Find periodic points of Tori %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Finding periodic points ... \n    ')
perPts = perPointFinder(tori, K, alpha);
toc
fprintf('\n')
plot(perPts(1,:),perPts(2,:),'.k')
%% Find rotation number %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Compute the rotation number... \n')
rhoMultiplier = 1;
diffRho = 1;
numPointsRho = 0;
rho = 0;
while (diffRho > 10^-rhoLimit && numPointsRho < 10^8) || rhoMultiplier < 25
    % Start with 100 points per tori and up the multiplier each iteration
    numPointsRho = K*1000*rhoMultiplier; 
    fprintf('Number of points per tori         %d\n', 1000*rhoMultiplier)
    trajectory2 = pointTrajectory(initialP, alpha, numPointsRho);
    tori2 = trajectorySeparator(trajectory2, K);
    rho(rhoMultiplier) = weightedBirkoffRotationNumber(tori2, K, perPts);
    fprintf('  Rho approximates as %.15g. \n',rho(rhoMultiplier))
    if rhoMultiplier > 1
        diffRho = abs(rho(rhoMultiplier) - rho(rhoMultiplier-1));
        fprintf('  Delta rho:          %e \n'...
            , abs(rho(rhoMultiplier)-rho(rhoMultiplier-1)))
    end
    fprintf('  Machine epsilon is  %e \n  ', eps(rho(rhoMultiplier)))
    fprintf(tableFile, ' %d & %.15f & %g & %g \\\\ \n', 1000*rhoMultiplier, rho(rhoMultiplier), diffRho, eps(rho(rhoMultiplier)));
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
fprintf(tableFile, '\\end{tabular} %% Matlab generated table');
diary off;
fclose all;