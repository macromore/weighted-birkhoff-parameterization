%% Tori hunter script

alpha = input('Alpha value? ');

figure
hold on
phasespacePlot(alpha);
trajectory = plot(0,0,'r.');
p = [1,0];

scale = 4;

axis([0, 2*pi, -scale, scale])

fprintf('Input (0,0) to stop.\n')

while p(1) ~= 0 || p(2) ~= 0
    p(1) = input('Input x-coordinate: ');
    p(2) = input('Input y-coordinate: ');
    fprintf('\n')
    if p(1) ~= 0 || p(2) ~= 0
        delete(trajectory);
        trajectory = trajectoryPlot(p, alpha, 10000, 'r');
    end
end