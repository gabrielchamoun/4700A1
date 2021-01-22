clc; close all;
set(0, 'DefaultFigureWindowStyle', 'docked')

eCount = 10;


tStep = 1e-12;
tStop = 0.01e-9;
time = 0;

eGroup = struct('x', 'y', 'vx', 'vy');
for i = 1 : eCount
    eGroup(i).x = rand()*200e-9;
    eGroup(i).y = rand()*100e-9;
    % Randomized velocities in direction and magnitude
    eGroup(i).vx = 1000*rand()*(2*randi([0 1])-1);
    eGroup(i).vy = 1000*rand()*(2*randi([0 1])-1);
end


counter = 1;
while time <= tStop

    for i = 1 : eCount
        plot(eGroup(i).x, eGroup(i).y)
        hold on
        
        % Updating Position
        eGroup(i).x(counter+1) = eGroup(i).x(counter) - eGroup(i).vx * tStep;
        eGroup(i).y(counter+1) = eGroup(i).y(counter) - eGroup(i).vy * tStep;

    end
    hold off
    time = time + tStep;
    counter = counter + 1;
    
end

