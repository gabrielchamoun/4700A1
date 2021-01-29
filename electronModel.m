clc; close all; clear all;
set(0, 'DefaultFigureWindowStyle', 'docked')


eCount = 1000;      % Total number of electrons
ePlotted = 10;      % Number of Electrons Plotted
tStop = 0.1e-9;     % Stop Time
dt = 1e-12;         % Time step 1ps

kB = 1.38066e-23;   % J/K
m0 = 4.66307e-26;   % Kg (Atomic mass of silicon = 28.0855)
mn = 0.26*m0;

Width = 200e-9;
Height = 100e-9;

% Thermal Velocity
Temp = 300; % K
vT = sqrt((kB*Temp)/mn);

eGroup = struct('x', 'y', 'vx', 'vy');  % Electron Object Classification
eVelocities = vT+rand(1,eCount);
eColours = rand(eCount,3);   % Random colours for plotting
for i = 1 : eCount
    eGroup(i).x = rand()*Width;
    eGroup(i).y = rand()*Height;
    % Randomized velocities in direction and magnitude ... 1000*rand()
    eGroup(i).vx = eVelocities(i)*(2*randi([0 1])-1);
    eGroup(i).vy = eVelocities(i)*(2*randi([0 1])-1);
end

plot(0,0);      % Init plot (comes up faster)
t = 0;          % Init time
counter = 1;    % Init Counter
while t < tStop

    for i = 1 : eCount
        
        % Updating Position
        eGroup(i).x(counter+1) = eGroup(i).x(counter) + eGroup(i).vx * dt;
        eGroup(i).y(counter+1) = eGroup(i).y(counter) + eGroup(i).vy * dt;
        
        %Plotting
        if (i <= 10)
            p = plot( [eGroup(i).x(counter), eGroup(i).x(counter+1)], ...
                [eGroup(i).y(counter), eGroup(i).y(counter+1)] );
            p.Color = eColours(i,:);
            hold on
        end
        % Conditional Statements
        % y = 200nm boundary
        if eGroup(i).x(counter+1) > Width
            eGroup(i).x(counter+1) = eGroup(i).x(counter+1) - Width;
        end
        
        % y = 0nm boundary
        if eGroup(i).x(counter+1) < 0
            eGroup(i).x(counter+1) = eGroup(i).x(counter+1) + Width;
        end
        
        % x = 100nm boundary
        if eGroup(i).y(counter+1) > Height
            diff = eGroup(i).y(counter+1) - Height;
            eGroup(i).y(counter+1) = Height - diff;
            eGroup(i).vy = -eGroup(i).vy;
        end
        
        % x = 0nm boundary
        if eGroup(i).y(counter+1) < 0
            diff = -eGroup(i).y(counter+1);
            eGroup(i).y(counter+1) = diff;
            eGroup(i).vy = -eGroup(i).vy;
        end
        
    end
    
    pause(0.001);    % To display trace as animation a delay is presented
    axis([0,Width,0,Height]);  % Plot Axis' being set
    t = t + dt;        % Incrementing Time
    counter = counter + 1;      % Incrementing Simulation counter
    
end
hold off
