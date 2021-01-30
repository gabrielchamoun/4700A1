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
Temp = 300;     % K
vT = sqrt((kB*Temp)/mn);

% scattering
% probability of scattering
tmin = 0.2e-10;
pScatter = 1 - exp(-dt/tmin);


% Initializing all electrons with a position and velocity
eGroup = struct('x', 'y', 'vx', 'vy');  % Electron Object Classification
eVelocities = vT + rand(1,eCount);
eColours = rand(eCount,3);   % Random colours for plotting
for i = 1 : eCount
    eGroup(i).x = rand()*Width;
    eGroup(i).y = rand()*Height;
    % Randomized velocities in direction and magnitude
    eGroup(i).vx = eVelocities(i)*(2*randi([0 1])-1);
    eGroup(i).vy = eVelocities(i)*(2*randi([0 1])-1);
end

t = 0;          % Init time
counter = 1;    % Init Counter
while t < tStop

    for i = 1 : eCount
        
        % Updating Position
        eGroup(i).x(counter+1) = eGroup(i).x(counter) + eGroup(i).vx * dt;
        eGroup(i).y(counter+1) = eGroup(i).y(counter) + eGroup(i).vy * dt;
        
        % scattering effect - Randomize direction/magnitude of velocity
        % probability of scattering based on p.
        if pScatter > rand()	% 'if true'
            eVelocities(i) = vT + rand();   % Regenerating a random velocity
            eGroup(i).vx = eVelocities(i)*(2*randi([0 1])-1);
            eGroup(i).vy = eVelocities(i)*(2*randi([0 1])-1);
        end
        
        
        % Plotting the first 10 electrons
        if (i <= 10)
            subplot(2,1,1)
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
    pause(0.001);       % Delay for animation
    axis([0,Width,0,Height]);  % Plot Axis' set
    t = t + dt;         % Incrementing Time
    
    Time(:,counter) = t;
    avgVelocity = mean(([eGroup(:).vx].^2 + [eGroup(:).vy].^2).^(1/2));
    Temp(:,counter) = ( (avgVelocity^2) * mn) / kB;
    subplot(2,1,2), plot(Time, Temp);
    
    counter = counter + 1;      % Incrementing Sim Counter
end
hold off
