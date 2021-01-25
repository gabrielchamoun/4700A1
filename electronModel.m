clc; close all;
set(0, 'DefaultFigureWindowStyle', 'docked')

eCount = 10;


tStep = 1e-12;
tStop = 0.1e-9;
time = 0;

eGroup = struct('x', 'y', 'vx', 'vy');
for i = 1 : eCount
    eGroup(i).x = rand()*200e-9;
    eGroup(i).y = rand()*100e-9;
    % Randomized velocities in direction and magnitude
    eGroup(i).vx = 1000*rand()*(2*randi([0 1])-1);
    eGroup(i).vy = 1000*rand()*(2*randi([0 1])-1);
end

plot(0,0);

counter = 1;
while time <= tStop

    for i = 1 : eCount
        %Plotting
        plot(eGroup(i).x, eGroup(i).y)
        hold on
        
        % Updating Position
        eGroup(i).x(counter+1) = eGroup(i).x(counter) - eGroup(i).vx * tStep;
        eGroup(i).y(counter+1) = eGroup(i).y(counter) - eGroup(i).vy * tStep;
        
        
        
        % Conditional Statements
        % y = 200nm boundary
        if eGroup(i).x(counter+1) > 200e-9
            eGroup(i).x(counter+1) = eGroup(i).x(counter+1) - 200e-9;
        end
        
        % y = 0nm boundary
        if eGroup(i).x(counter+1) < 0
            eGroup(i).x(counter+1) = eGroup(i).x(counter+1) + 200e-9;
        end
        
        % x = 100nm boundary
        if eGroup(i).y(counter+1) > 100e-9
            diff = eGroup(i).y(counter+1) - 100e-9;
            eGroup(i).y(counter+1) = 100e-9 - diff;
            eGroup(i).vy = -eGroup(i).vy;
        end
        
        % x = 0nm boundary
        if eGroup(i).y(counter+1) < 0
            diff = -eGroup(i).y(counter+1);
            eGroup(i).y(counter+1) = diff;
            eGroup(i).vy = -eGroup(i).vy;
        end
        
    end
    hold off
    pause(0.001)    % To display trace as animation a delay is presented
    axis([0,200e-9,0,100e-9]);  % Plot Axis' being set
    time = time + tStep;        % Incrementing Time
    counter = counter + 1;      % Incrementing Simulation counter
    
end

